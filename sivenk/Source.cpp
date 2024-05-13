#include "Base.h"
#include "Generate.h"

int main()
{
	Make_grid("");
	Create_time_grid();
	FEM task;
	std::vector<double> res;
	task.SolveTask_time();
	//task.SolveTask(res);
	return 0;
}
