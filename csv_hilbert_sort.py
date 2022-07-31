#!/usr/bin/env python3
'''
Simple command line tool for reading in a csv that contains columns named
'x' and 'y', calculating their Hilbert index, and sorting the entire file by the
Hilbert index.
'''

import argparse
import sys
import pandas as pd
from hilbertcurve.hilbertcurve import HilbertCurve


def calc_scaled_xy(df, scale_factor=1):
    '''
    Calculate the scaled (to 0...scale_factor) x and y coordinates for the dataframe.
    '''
    xmin, ymin = 0, -90
    xmax, ymax = 360, 90

    x_range = xmax - xmin
    y_range = ymax - ymin

    df['scaled_x'] = scale_factor * (df['ra'] - xmin)/x_range
    df['scaled_y'] = scale_factor * (df['dec'] - ymin)/y_range
    return df

def calc_hilbert_distance(hilbert_curve, df, scaled_coords):
    '''
    calculate the distance along the hilbert curve for each point in the dataframe.
    '''
    floatpoints = df[scaled_coords].values.tolist()
    intpoints = [[round(item[0]), round(item[1])] for item in floatpoints]
    df['hilbert_index'] = hilbert_curve.distances_from_points(intpoints)
    return df

def hilbert_sort(df, hilbert_order=10, rescale=True):
    '''
    Sort the dataframe by the Hilbert index.
    '''
    hc = HilbertCurve(hilbert_order, 2)
    scale_factor = 2**hilbert_order - 1

    if rescale:
        df = calc_scaled_xy(df, scale_factor)
        scaled_coords = ['scaled_x', 'scaled_y']
    else:
        scaled_coords = ['x', 'y']
    df = calc_hilbert_distance(hc, df, scaled_coords)
    df = df.sort_values(by='hilbert_index')

    return df

def csv_hilbert_sort():
    parser = argparse.ArgumentParser(description="Sort (x,y) pairs by position along hilbert curve.")
    parser.add_argument('input_csv', help="Input csv file (required)")
    parser.add_argument('--output', '-o', type=argparse.FileType('wr'), default=sys.stdout,
                        help="Output csv file (default: stdout)")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="Maintain calculated columns in output csv")
    parser.add_argument('--keep_scale', '-k', action='store_true',
                        help="Don't rescale values")
    parser.add_argument('--hilbert_order', '-n', default=10, type=int,
                        help="Order of the hilbert curve (default: 10)")

    args = parser.parse_args()
    df = pd.read_csv(args.input_csv, low_memory=False)
    df.insert(0, 'hilbert_index', 0)
    df = hilbert_sort(df, args.hilbert_order, not args.keep_scale)
    if not args.verbose:
        df.drop(['hilbert_index'], axis=1, inplace=True)
    if not args.keep_scale:
        df.drop(['scaled_x', 'scaled_y'], axis=1, inplace=True)
    pd.DataFrame(df).to_csv(args.output, index=False)


if __name__ == "__main__":
    csv_hilbert_sort()
