# -*- coding: utf-8 -*-
from model_g10_v2 import *
from run_manager import *

"""
This script tests in CPU time how long does one 
full model take: 1 to 355 days.

It has a lot of boiler plate code not relevant here 
because it is just a copy of the "coupling" tests.
This should not really affect the computation time, 
such that such additional lines can be expected to 
be of order n, T(n) = O(n).
"""


def test():
    origin = os.getcwd()
    print(origin)
    """
    Input parameters to start BEACH
    """
    names = ['z3_factor',
             'cZ0Z1',
             'cZ',
             'c_adr',
             'k_g',
             'gamma01', 'gammaZ',
             'f_transp',
             'f_evap',
             'f_oc', 'k_oc',
             'beta_runoff',
             'dt_50_aged',
             'dt_50_ab',
             'dt_50_ref',
             'epsilon_iso',
             'beta_moisture'
             ]
    best_values = np.loadtxt("best_vector.txt")  # Return a vector, with same values as names
    upper = np.ones(len(best_values)).tolist()  # Just a vector of 1's

    # Reduce "best" drainge coefficient
    correct_factor = .5*.5*.5
    new_values = []
    for i in range(len(names)):
        if i != names.index("c_adr"):
            new_values.append(best_values[i])
        else:
            new_values.append(best_values[i]*correct_factor)

    # First BEACH run
    couple = False  # If False, run a pure BEACH run, without coupling.

    firstTimeStep = 1  # 166 -> 14/03/2016
    runs = 0
    correct = False
    for period in range(len(lisem_runs)):
        os.chdir(origin)
        if runs == 0:
            runs += 1
            firstTimeStep = firstTimeStep
            if couple:
                lastTimeStep = lisem_runs[period]
                param_values = new_values
            else:
                lastTimeStep = 355
                param_values = best_values
        else:
            runs += 1
            firstTimeStep = lisem_runs[period - 1] - 1
            lastTimeStep = lisem_runs[period]

        print("Starting BEACH")
        print("run = " + str(runs))
        print("period = " + str(period))
        myAlteck16 = BeachModel("maps\\static\\clone_nom.map",
                                names, param_values, upper, period,
                                lisem_runs, end_days_remain, event_type,
                                map_list, hydro_row, ts_new,
                                staticDT50=False, correction=correct, coupled=couple)
        dynamicModel = DynamicFramework(myAlteck16,
                                        firstTimestep=firstTimeStep,
                                        lastTimeStep=lastTimeStep  # -> 183
                                        )  # an instance of the Dynamic Framework
        dynamicModel.run()
        end_days_remain.pop(0)
        correct = True  # Will tell BEACH to import LISEM outputs instead of initial conditions on next run
        print("Finished")
        print(datetime.today().strftime('%Y-%m-%d %HH:%MM'))
        if not couple:
            break
        else:
            print("Preparing LISEM")
            # Modify run file (appropriate time-steps for recording maps)
            # total: 510 minutes; timestep: 1248 sec/timestep ,
            #   break-day: 208min x 60sec = 12480 sec;
            #       output interval = 12480/1248 = 10.

            os.chdir(os.getcwd() + r"\\LISEM")
            path_run = r"..\\Alteck" + str(lisem_runs[period]) + ".run"
            command_line = r"lisem -b -r " + str(path_run)
            args = shlex.split(command_line)

            print("Starting LISEM")
            code = subprocess.call(args)

            if code == -1073741819:
                # end_run = False
                print("Code was good: " + str(code))
                if runs == len(lisem_runs):  # Do final run!
                    firstTimeStep = lisem_runs[period] - 1
                    # end_run = True
                    runs += 1
                    print("Starting last BEACH run")
                    print("run = " + str(runs))
                    print("period = " + str(period+1))
                    os.chdir(origin)
                    myAlteck16 = BeachModel("maps\\static\\clone_nom.map",
                                            names, param_values, upper, period + 1,
                                            lisem_runs, end_days_remain, event_type,
                                            map_list, hydro_row, ts_new,
                                            staticDT50=False, correction=correct, coupled=couple)
                    dynamicModel = DynamicFramework(myAlteck16,
                                                    firstTimestep=firstTimeStep,
                                                    lastTimeStep=286  # -> 183
                                                    )  # an instance of the Dynamic Framework
                    dynamicModel.run()
                    # if len(end_days_remain) > 0:
                    # end_days_remain.pop(0)
                    break
                else:
                    continue
            else:
                print("Error code: " + str(code))
                print("Exiting coupled loop")
                break

    print("All runs finished")
    print("Remaining days to couple: ", end_days_remain)


if __name__ == '__main__':
    import timeit
    print(timeit.timeit("test()", setup="from __main__ import test", number=1))
