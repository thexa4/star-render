#include "star_record.hpp"
#include <iostream>
#include <queue>
#include <climits>
#include <cerrno>
#include <sys/resource.h>

class MergeState {
	FILE *file;

public:
	MergeState(const char *filename) {
		done = false;
		failed = false;
		file = fopen(filename, "r");
		if (file == NULL) {
			std::cerr << "Error while opening " << filename << ": " << errno << std::endl;
			done = true;
			failed = true;
		} else {
			ReadSample();
			this->filename = filename;
		}
	}

	Star current_sample;
	bool done;
	bool failed;
	const char *filename;

	void ReadSample() {
		if (done)
			return;

		if (fread(&current_sample, sizeof(current_sample), 1, file) == 0) {
			fclose(file);
			done = true;
		}
	}
};

bool operator < (const MergeState& lhs, const MergeState& rhs)
{
    return lhs.current_sample.hilbert_index > rhs.current_sample.hilbert_index;
}

bool operator > (const MergeState& lhs, const MergeState& rhs)
{
    return lhs.current_sample.hilbert_index < rhs.current_sample.hilbert_index;
}

int main(int argc, char *argv[]) {

	rlimit lim;
	getrlimit(RLIMIT_NOFILE, &lim);
	lim.rlim_cur = std::max((int)lim.rlim_max, argc + 10);
	setrlimit(RLIMIT_NOFILE, &lim);

	std::priority_queue<MergeState> states;
	for (int i = 1; i < argc; i++) {
		MergeState s{argv[i]};
		if (s.failed)
			return errno;
		states.push(s);
	}

	while(!states.empty()) {

		MergeState s = states.top();
		states.pop();
		auto next = states.top();

		while (!s.done && s.current_sample.hilbert_index <= next.current_sample.hilbert_index) {
			fwrite(&s.current_sample, sizeof(s.current_sample), 1, stdout);
			s.ReadSample();
		}

		if (!s.done)
			states.push(s);
	}

	return 0;
}
