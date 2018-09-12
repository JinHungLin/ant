#include <iostream>
#include <random>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "beta.h"

using namespace std;

#define CT 10 // maximum time for a work station
#define FMPP 0.5 // follow the max pheromone path probability

enum KItem
{
	BloodCheck,
	AbdUltra,
	Endoscope,
	NumKItem
};

enum PItem
{
	EyeCheck,
	HeartCheck,
	XRay,
	BoneDensity,
	BasicCheck,
	HeartUltra,
	NeckUltra,
	ThyroidUltra,
	ProstateUltra,
	GynecologyUltra,
	BreastUltra,
	MRI,
	CTScan,
	LungTumorScreening,
	CoronaryCalcificationAssessment,
	PeripheralStenosis,
	PulmonaryCheck,
	Radiology,
	Ophthalmologist,
	FamilyMedicine,
	NumPItem
};

enum ItemType
{
	K_TYPE,
	P_TYPE
};

enum FindTaskStatus
{
	TASK_LIST_EMPTY = -2,
    WORK_STATION_OVER_TIME = -1
};

struct Task
{
	int type;
	int index;
	bool need_schedule;
	bool schedule_done;
	double task_time;
	struct Task *next_task;
};

struct Model
{
	int people;
	int num_of_tasks;
	struct Task **task_list;
};

struct WorkStation
{
	int num_of_tasks;
	struct Task *first_task;
	double work_station_time;
	struct WorkStation *next_work_station;
};

struct Solution
{
    int num_of_work_stations;
	struct WorkStation *first_work_station;
	double objective_function;
};

double triangular(double a, double b, double c) {
    //srand (time(NULL));
    double U = rand() / (double) RAND_MAX;
    double F = (c - a) / (b - a);
    if (U <= F)
        return a + sqrt(U * (b - a) * (c - a));
    else
        return b - sqrt((1 - U) * (b - a) * (b - c));
}

double Pheromone_KK [NumKItem][NumKItem];
double Pheromone_KP [NumKItem][NumPItem];
double Pheromone_PP [NumPItem][NumPItem];

struct Task K_task [NumKItem];
struct Task P_task [NumPItem];

struct Solution* build_feasible_solution(void);
int find_K_task(struct Task *currentTask, double work_station_time);
double calulate_task_time(int type, int index);

int main()
{
    srand (time(NULL));

	int N = 50;
	int SC = 200;
	double a = 0.2;
	double b0 = 8.0;

	// initialize pheromone
	for (int i=0; i<NumKItem; i++)
        for (int j=0; j<NumKItem; j++)
            Pheromone_KK[i][j] = b0;

    for (int i=0; i<NumKItem; i++)
        for (int j=0; j<NumPItem; j++)
            Pheromone_KP[i][j] = b0;

    for (int i=0; i<NumPItem; i++)
        for (int j=0; j<NumPItem; j++)
            Pheromone_PP[i][j] = b0;



	for (int i=0; i<NumKItem; i++)
    {
        K_task[i].type = K_TYPE;
        K_task[i].index = i;
        K_task[i].need_schedule = 0;
        K_task[i].schedule_done = 0;
        K_task[i].task_time = 0;
        K_task[i].next_task = nullptr;
    }
	for (int i=0; i<NumPItem; i++)
    {
        P_task[i].type = P_TYPE;
        P_task[i].index = i;
        P_task[i].need_schedule = 0;
        P_task[i].schedule_done = 0;
        P_task[i].task_time = 0;
        P_task[i].next_task = nullptr;
    }

    struct Model model[10];

    model[0].num_of_tasks = 7;
    model[1].num_of_tasks = 8;
    model[2].num_of_tasks = 12;
    model[3].num_of_tasks = 10;
    model[4].num_of_tasks = 9;
    model[5].num_of_tasks = 9;
    model[6].num_of_tasks = 9;
    model[7].num_of_tasks = 9;
    model[8].num_of_tasks = 20;
    model[9].num_of_tasks = 20;

    for (int i=0; i<10; i++)
    {
        struct Task **tk = new struct Task* [model[i].num_of_tasks];
        model[i].task_list = tk;
    }

    model[0].task_list[0] = &K_task[BloodCheck];
    model[0].task_list[1] = &K_task[AbdUltra];
    model[0].task_list[2] = &P_task[EyeCheck];
    model[0].task_list[3] = &P_task[HeartCheck];
    model[0].task_list[4] = &P_task[XRay];
    model[0].task_list[5] = &P_task[BasicCheck];
    model[0].task_list[6] = &P_task[FamilyMedicine];
    model[1].task_list[0] = &K_task[BloodCheck];
    model[1].task_list[1] = &K_task[AbdUltra];
    model[1].task_list[2] = &P_task[EyeCheck];
    model[1].task_list[3] = &P_task[HeartCheck];
    model[1].task_list[4] = &P_task[XRay];
    model[1].task_list[5] = &P_task[BoneDensity];
    model[1].task_list[6] = &P_task[BasicCheck];
    model[1].task_list[7] = &P_task[FamilyMedicine];
    model[2].task_list[0] = &K_task[BloodCheck];
    model[2].task_list[1] = &K_task[AbdUltra];
    model[2].task_list[2] = &P_task[EyeCheck];
    model[2].task_list[3] = &P_task[HeartCheck];
    model[2].task_list[4] = &P_task[XRay];
    model[2].task_list[5] = &P_task[BoneDensity];
    model[2].task_list[6] = &P_task[BasicCheck];
    model[2].task_list[7] = &P_task[HeartUltra];
    model[2].task_list[8] = &P_task[NeckUltra];
    model[2].task_list[9] = &P_task[MRI];
    model[2].task_list[10] = &P_task[CoronaryCalcificationAssessment];
    model[2].task_list[11] = &P_task[FamilyMedicine];
    model[3].task_list[0] = &K_task[BloodCheck];
    model[3].task_list[1] = &K_task[AbdUltra];
    model[3].task_list[2] = &P_task[EyeCheck];
    model[3].task_list[3] = &P_task[HeartCheck];
    model[3].task_list[4] = &P_task[XRay];
    model[3].task_list[5] = &P_task[BoneDensity];
    model[3].task_list[6] = &P_task[BasicCheck];
    model[3].task_list[7] = &P_task[NeckUltra];
    model[3].task_list[8] = &P_task[MRI];
    model[3].task_list[9] = &P_task[FamilyMedicine];
    model[4].task_list[0] = &K_task[BloodCheck];
    model[4].task_list[1] = &K_task[AbdUltra];
    model[4].task_list[2] = &K_task[Endoscope];
    model[4].task_list[3] = &P_task[EyeCheck];
    model[4].task_list[4] = &P_task[HeartCheck];
    model[4].task_list[5] = &P_task[XRay];
    model[4].task_list[6] = &P_task[BasicCheck];
    model[4].task_list[7] = &P_task[MRI];
    model[4].task_list[8] = &P_task[FamilyMedicine];
    model[5].task_list[0] = &K_task[BloodCheck];
    model[5].task_list[1] = &K_task[AbdUltra];
    model[5].task_list[2] = &P_task[EyeCheck];
    model[5].task_list[3] = &P_task[HeartCheck];
    model[5].task_list[4] = &P_task[XRay];
    model[5].task_list[5] = &P_task[BasicCheck];
    model[5].task_list[6] = &P_task[MRI];
    model[5].task_list[7] = &P_task[LungTumorScreening];
    model[5].task_list[8] = &P_task[FamilyMedicine];
    model[6].task_list[0] = &K_task[BloodCheck];
    model[6].task_list[1] = &K_task[AbdUltra];
    model[6].task_list[2] = &P_task[EyeCheck];
    model[6].task_list[3] = &P_task[HeartCheck];
    model[6].task_list[4] = &P_task[XRay];
    model[6].task_list[5] = &P_task[BoneDensity];
    model[6].task_list[6] = &P_task[BasicCheck];
    model[6].task_list[7] = &P_task[LungTumorScreening];
    model[6].task_list[8] = &P_task[FamilyMedicine];
    model[7].task_list[0] = &K_task[BloodCheck];
    model[7].task_list[1] = &K_task[AbdUltra];
    model[7].task_list[2] = &K_task[Endoscope];
    model[7].task_list[3] = &P_task[EyeCheck];
    model[7].task_list[4] = &P_task[HeartCheck];
    model[7].task_list[5] = &P_task[XRay];
    model[7].task_list[6] = &P_task[BoneDensity];
    model[7].task_list[7] = &P_task[BasicCheck];
    model[7].task_list[8] = &P_task[FamilyMedicine];
    model[8].task_list[0] = &K_task[BloodCheck];
    model[8].task_list[1] = &K_task[AbdUltra];
    model[8].task_list[2] = &P_task[EyeCheck];
    model[8].task_list[3] = &P_task[HeartCheck];
    model[8].task_list[4] = &P_task[XRay];
    model[8].task_list[5] = &P_task[HeartUltra];
    model[8].task_list[6] = &P_task[BasicCheck];
    model[8].task_list[7] = &P_task[NeckUltra];
    model[8].task_list[8] = &P_task[ThyroidUltra];
    model[8].task_list[9] = &P_task[ProstateUltra];
    model[8].task_list[10] = &P_task[GynecologyUltra];
    model[8].task_list[11] = &P_task[BreastUltra];
    model[8].task_list[12] = &P_task[MRI];
    model[8].task_list[13] = &P_task[CTScan];
    model[8].task_list[14] = &P_task[LungTumorScreening];
    model[8].task_list[15] = &P_task[PeripheralStenosis];
    model[8].task_list[16] = &P_task[PulmonaryCheck];
    model[8].task_list[17] = &P_task[Radiology];
    model[8].task_list[18] = &P_task[Ophthalmologist];
    model[8].task_list[19] = &P_task[FamilyMedicine];
    model[9].task_list[0] = &K_task[BloodCheck];
    model[9].task_list[1] = &K_task[AbdUltra];
    model[9].task_list[2] = &P_task[EyeCheck];
    model[9].task_list[3] = &P_task[HeartCheck];
    model[9].task_list[4] = &P_task[XRay];
    model[9].task_list[5] = &P_task[HeartUltra];
    model[9].task_list[6] = &P_task[BasicCheck];
    model[9].task_list[7] = &P_task[NeckUltra];
    model[9].task_list[8] = &P_task[ThyroidUltra];
    model[9].task_list[9] = &P_task[ProstateUltra];
    model[9].task_list[10] = &P_task[GynecologyUltra];
    model[9].task_list[11] = &P_task[BreastUltra];
    model[9].task_list[12] = &P_task[MRI];
    model[9].task_list[13] = &P_task[CoronaryCalcificationAssessment];
    model[9].task_list[14] = &P_task[LungTumorScreening];
    model[9].task_list[15] = &P_task[PeripheralStenosis];
    model[9].task_list[16] = &P_task[PulmonaryCheck];
    model[9].task_list[17] = &P_task[Radiology];
    model[9].task_list[18] = &P_task[Ophthalmologist];
    model[9].task_list[19] = &P_task[FamilyMedicine];

	for (int i=0; i<10; i++)
    {
        cout << "Please enter the # of people in model " << i << ": ";
        cin >> model[i].people;
        if (model[i].people != 0)
            for (int j=0; j<model[i].num_of_tasks; j++)
                model[i].task_list[j]->need_schedule = 1;
    }

	return 0;
}

struct Solution* build_feasible_solution(void)
{
    bool K_type_task_empty = 0;
    bool P_type_task_empty = 0;

    for (int i=0; i<NumKItem; i++)
        K_task[i].task_time = calulate_task_time(K_TYPE, i);

    for (int i=0; i<NumPItem; i++)
        P_task[i].task_time = calulate_task_time(P_TYPE, i);

	// create a solution
	struct Solution* sol = new struct Solution;

	// open work station
    struct WorkStation* ws = new struct WorkStation;
    *ws = {
        .num_of_tasks = 0,
        .first_task = nullptr,
        .work_station_time = 0,
        .next_work_station = nullptr
    };

    sol->first_work_station = ws;

    struct Task *tk = nullptr;

    // arrange K tasks
    while (true)
    {
        int tk_index = find_K_task(tk, ws->work_station_time);
        if (tk_index == TASK_LIST_EMPTY)
        {
            break;
        }
        else if (tk_index == WORK_STATION_OVER_TIME)
        {
            // create new work station

            continue;
        }
        else
        {
            // arrange the K task found
        }
	}
}

int find_K_task(struct Task *currentTask, double work_station_time)
{
	int total_task_need_schedule = 0;
	int total_task_over_time = 0;

    for (int i=0; i<NumKItem; i++)
	{
        if (K_task[i].need_schedule == 1 && K_task[i].schedule_done != 1)
            total_task_need_schedule++;
		if (K_task[i].need_schedule == 1 && K_task[i].schedule_done != 1 && work_station_time + K_task[i].task_time > CT)
			total_task_over_time++;
	}

    if (total_task_need_schedule == 0)
        return TASK_LIST_EMPTY;

	if (total_task_need_schedule == total_task_over_time)
		return WORK_STATION_OVER_TIME;

	double p = rand() / (double) RAND_MAX;
    if (currentTask == nullptr || p > FMPP)
    {
        // random selection
        double sel_prob = 1.0 / (double) (total_task_need_schedule - total_task_over_time);
        for (int i=0; i<NumKItem; i++)
        {
            if (K_task[i].need_schedule == 1 && K_task[i].schedule_done != 1 && work_station_time + K_task[i].task_time <= CT)
            {
				p = rand() / (double) RAND_MAX;
				if (p < sel_prob)
				{
					currentTask->next_task = &K_task[i];
					K_task[i].schedule_done = 1;
					return i;
				}
            }
        }
    }
    else
    {
        // select max pheromone
        double max_Pheromone = 0;
        int max_Pheromone_Task_Index;
        for (int i=0; i<NumKItem; i++)
        {
            if (K_task[i].need_schedule == 1 && K_task[i].schedule_done != 1 && work_station_time + K_task[i].task_time <= CT)
            {
                if (Pheromone_KK[currentTask->index][i] > max_Pheromone)
                {
                    max_Pheromone = Pheromone_KK[currentTask->index][i];
                    max_Pheromone_Task_Index = i;
                }
            }
        }
        currentTask->next_task = &K_task[max_Pheromone_Task_Index];
        K_task[max_Pheromone_Task_Index].schedule_done = 1;
        return max_Pheromone_Task_Index;
    }
}

int find_P_task(struct Task *currentTask, double work_station_time)
{
	int total_task_need_schedule = 0;
	int total_task_over_time = 0;

    for (int i=0; i<NumPItem; i++)
	{
        if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1)
            total_task_need_schedule++;
		if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1 && work_station_time + P_task[i].task_time > CT)
			total_task_over_time++;
	}

    if (total_task_need_schedule == 0)
        return TASK_LIST_EMPTY;

	if (total_task_need_schedule == total_task_over_time)
		return WORK_STATION_OVER_TIME;

	double p = rand() / (double) RAND_MAX;
    if (currentTask == &K_task[NumKItem] && p > FMPP)
    {
        // random selection
        double sel_prob = 1.0 / (double) (total_task_need_schedule - total_task_over_time);
        for (int i=0; i<NumPItem; i++)
        {
            if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1 && work_station_time + P_task[i].task_time <= CT)
            {
				p = rand() / (double) RAND_MAX;
				if (p < sel_prob)
				{
					currentTask->next_task = &P_task[i];
					P_task[i].schedule_done = 1;
					return i;
				}
            }
        }
    }
    else if (currentTask == &K_task[NumKItem] && p <= FMPP)
    {
        // select max pheromone
        double max_Pheromone = 0;
        int max_Pheromone_Task_Index;
        for (int i=0; i<NumKItem; i++)
        {
            if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1 && work_station_time + P_task[i].task_time <= CT)
            {
                if (Pheromone_KP[currentTask->index][i] > max_Pheromone)
                {
                    max_Pheromone = Pheromone_KP[currentTask->index][i];
                    max_Pheromone_Task_Index = i;
                }
            }
        }
        currentTask->next_task = &P_task[max_Pheromone_Task_Index];
        P_task[max_Pheromone_Task_Index].schedule_done = 1;
        return max_Pheromone_Task_Index;
    }
    else if (currentTask == &P_task[NumPItem] && p > FMPP)
    {
        // random selection
        double sel_prob = 1.0 / (double) (total_task_need_schedule - total_task_over_time);
        for (int i=0; i<NumPItem; i++)
        {
            if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1 && work_station_time + P_task[i].task_time <= CT)
            {
				p = rand() / (double) RAND_MAX;
				if (p < sel_prob)
				{
					currentTask->next_task = &P_task[i];
					P_task[i].schedule_done = 1;
					return i;
				}
            }
        }
    }
    else if (currentTask == &P_task[NumPItem] && p <= FMPP)
    {
        // select max pheromone
        double max_Pheromone = 0;
        int max_Pheromone_Task_Index;
        for (int i=0; i<NumKItem; i++)
        {
            if (P_task[i].need_schedule == 1 && P_task[i].schedule_done != 1 && work_station_time + P_task[i].task_time <= CT)
            {
                if (Pheromone_PP[currentTask->index][i] > max_Pheromone)
                {
                    max_Pheromone = Pheromone_PP[currentTask->index][i];
                    max_Pheromone_Task_Index = i;
                }
            }
        }
        currentTask->next_task = &P_task[max_Pheromone_Task_Index];
        P_task[max_Pheromone_Task_Index].schedule_done = 1;
        return max_Pheromone_Task_Index;
    }
}

double calulate_task_time(int type, int index)
{
	double task_time;
	default_random_engine e;

	if (type == K_TYPE)
	{
		switch(index)
		{
			case BloodCheck:
            {
				gamma_distribution<double> distribution(0.438, 3);
				task_time = 1.5 + distribution(e);
				break;
            }
			case AbdUltra:
            {
				lognormal_distribution<double> distribution(3.11, 2.25);
                task_time = 4.5 + distribution(e);
				break;
            }
			case Endoscope:
            {
				sftrabbit::beta_distribution<> distribution(0.373, 0.373);
				task_time = 15.5 + 11 * distribution(e);
				break;
            }
		}
	}
	else if (type == P_TYPE)
	{
		switch(index)
		{
			case EyeCheck:
            {
                task_time = triangular(2.5, 3.5, 4.5);
				break;
            }
			case HeartCheck:
            {
				sftrabbit::beta_distribution<> distribution(2.21, 2.27);
				task_time = 3.5 + 3 * distribution(e);
				break;
            }
			case XRay:
            {
				weibull_distribution<double> distribution(1.39, 3.72);
                task_time = 0.5 + distribution(e);
				break;
            }
			case BoneDensity:
            {
				uniform_real_distribution<double> distribution(1.5, 4.5);
                task_time = distribution(e);
				break;
            }
			case BasicCheck:
            {
				weibull_distribution<double> distribution(0.854, 1.9);
                task_time = 2.5 + distribution(e);
				break;
            }
			case HeartUltra:
            {
				sftrabbit::beta_distribution<> distribution(0.0607, 0.0819);
				task_time = 11.5 + 9 * distribution(e);
				break;
            }
			case NeckUltra:
            {
				sftrabbit::beta_distribution<> distribution(0.0804, 0.0982);
				task_time = 13.5 + 10 * distribution(e);
				break;
            }
			case ThyroidUltra:
            {
				sftrabbit::beta_distribution<> distribution(0.741, 1.07);
				task_time = 4.5 + 10 * distribution(e);
				break;
            }
			case ProstateUltra:
            {
				sftrabbit::beta_distribution<> distribution(0.741, 1.07);
				task_time = 4.5 + 10 * distribution(e);
				break;
            }
			case GynecologyUltra:
            {
				sftrabbit::beta_distribution<> distribution(1.08, 1.08);
				task_time = 2.5 + 6 * distribution(e);
				break;
            }
			case BreastUltra:
            {
				sftrabbit::beta_distribution<> distribution(1.08, 1.08);
				task_time = 2.5 + 6 * distribution(e);
				break;
            }
			case MRI:
            {
				sftrabbit::beta_distribution<> distribution(0.654, 0.609);
				task_time = 16.5 + 9 * distribution(e);
				break;
            }
			case CTScan:
            {
				sftrabbit::beta_distribution<> distribution(0.654, 0.609);
				task_time = 16.5 + 9 * distribution(e);
				break;
            }
			case LungTumorScreening:
            {
				exponential_distribution<double> distribution(1.3);
                task_time = 2.5 + distribution(e);
				break;
            }
			case CoronaryCalcificationAssessment:
            {
				exponential_distribution<double> distribution(1.3);
                task_time = 2.5 + distribution(e);
				break;
            }
			case PeripheralStenosis:
            {
				weibull_distribution<double> distribution(0.681, 2.11);
                task_time = 4.5 + distribution(e);
				break;
            }
			case PulmonaryCheck:
            {
				sftrabbit::beta_distribution<> distribution(0.449, 0.508);
				task_time = 8.5 + 16 * distribution(e);
				break;
            }
			case Radiology:
            {
                task_time = triangular(5, 10, 15);
				break;
            }
			case Ophthalmologist:
            {
				sftrabbit::beta_distribution<> distribution(0.794, 0.961);
				task_time = 1.5 + 3 * distribution(e);
				break;
            }
			case FamilyMedicine:
            {
				lognormal_distribution<double> distribution(3.47, 2.99);
                task_time = 3.5 + distribution(e);
				break;
            }
		}
	}
}
