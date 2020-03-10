import traceback
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.grids import add_flo2d_initial_conditions, get_flo2d_initial_conditions, \
    clear_initial_conditions
from db_adapter.constants import CURW_SIM_HOST, CURW_SIM_PORT, CURW_SIM_USERNAME, CURW_SIM_PASSWORD, CURW_SIM_DATABASE
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150

flo2d_model_list = [FLO2D_250, FLO2D_150]

print(" Add flo2d initial conditions")

try:

    pool = get_Pool(host=CURW_SIM_HOST, port=CURW_SIM_PORT, user=CURW_SIM_USERNAME, password=CURW_SIM_PASSWORD, db=CURW_SIM_DATABASE)

    for flo2d_model in flo2d_model_list:
        print("Clear existing {} initial conditions".format(flo2d_model))
        print("{} initial conditions cleared.".format(clear_initial_conditions(pool=pool, flo2d_model=flo2d_model)))

        print("Add {} initial conditions".format(flo2d_model))
        add_flo2d_initial_conditions(pool=pool, flo2d_model=flo2d_model,
                                     initial_condition_file_path='grid_maps/flo2d/initial_conditions/{}_initial_cond.csv'
                                     .format(flo2d_model))
        print("{} initial conditions added".format(len(get_flo2d_initial_conditions(pool=pool, flo2d_model=flo2d_model))))

except Exception as e:
    traceback.print_exc()
finally:
    destroy_Pool(pool=pool)
    print("Process Finished.")
