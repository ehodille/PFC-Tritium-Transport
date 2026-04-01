from scenario import Scenario, Pulse
import pandas as pd
from plasma_data_handling import PlasmaDataHandling

data_folder = "data"
plasma_data_handling = PlasmaDataHandling(
    pulse_type_to_data={
        "FP": pd.read_csv(data_folder + "/Binned_Flux_Data_new.dat", delimiter=","),
        "FP_D": pd.read_csv(data_folder + "/Binned_Flux_Data_just_D_pulse.dat", delimiter=",", comment='#'),
        "ICWC": pd.read_csv(data_folder + "/ICWC_data.dat", delimiter=","),
        "GDC": pd.read_csv(data_folder + "/GDC_data.dat", delimiter=","),
    },
    path_to_RISP_data=data_folder + "/RISP_data",
    path_to_ROSP_data=data_folder + "/ROSP_data",
    path_to_RISP_wall_data=data_folder + "/RISP_Wall_data.dat",
)

# Example pulse
example_fp = Pulse(
    pulse_type="FP",
    nb_pulses=1,
    ramp_up=1000,
    steady_state=2000,
    ramp_down=1000,
    waiting=2000,
    tritium_fraction=0.2,
    heat_scaling=0.33,
    flux_scaling=0.25,
)

scenario = Scenario(pulses=[
    example_fp,
])
