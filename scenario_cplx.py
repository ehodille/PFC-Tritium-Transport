import pandas as pd
from typing import List
import warnings

class Pulse:
    pulse_type: str
    nb_pulses: int
    steady_state: float
    steady_STATE: List[float]
    transition: List[float]
    fraction_ss: List[float]
    ramp_up: float
    ramp_down: float
    waiting: float
    tritium_fraction: float # tritium fraction = T/(D+T)
    heat_scaling: float
    flux_scaling: float

    def __init__(
        self,
        pulse_type: str,
        nb_pulses: int,
        transition:List[float],
        steady_STATE:List[float],
        fraction_ss:List[float],
        waiting: float,
        tritium_fraction: float,  # tritium fraction = T/(D+T)
        heat_scaling: float=1.0,  # scaling factor for heat loads
        flux_scaling: float=1.0,  # scaling factor for particle fluxes
    ):
        self.pulse_type = pulse_type
        self.nb_pulses = nb_pulses
        self.steady_STATE = steady_STATE
        if len(transition) == len(steady_STATE)+1:
            self.transition = transition
            self.ramp_up = transition[0]        # first value of transition
            self.ramp_down = transition[1]      # last value of transition
            self.steady_state = sum(steady_STATE)+sum(transition[1:-1])
        else:
            print('ERROR: len(transition) != len(steady_STATE)+1')
            sys.exit('exiting ...')
        if len(fraction_ss) == len(steady_STATE):
            self.fraction_ss = fraction_ss
        self.waiting = waiting
        self.tritium_fraction = tritium_fraction
        self.heat_scaling = heat_scaling
        self.flux_scaling = flux_scaling

    @property
    def total_duration(self) -> float: 

        all_zeros = (
            sum(self.transition) == 0
            and sum(self.steady_STATE) == 0
            and self.waiting == 0
        )
        if self.pulse_type == "RISP" and all_zeros:
            msg = "RISP pulse has all zeros for ramp_up, steady_state, ramp_down, waiting. "
            msg += "Setting hardcoded values. Please check the values in the scenario file."
            warnings.warn(msg, UserWarning)
            
            self.transition = [10,10]
            self.ramp_up = 10
            self.steady_state = 250
            self.steady_STATE = [250]
            self.waiting = 1530

        return sum(self.transition) + sum(self.steady_STATE) + self.waiting

    @property
    def duration_no_waiting(self) -> float:
        return self.total_duration - self.waiting


class Scenario:
    def __init__(self, pulses: List[Pulse] = None, baking_temp: float = None):
        """Initializes a Scenario object containing several pulses.

        Args:
            pulses: The list of pulses in the scenario. Each pulse is a Pulse object.
            baking_temp: Temperature (K) during baking pulses. Required if any
                pulse has pulse_type="BAKE".
        """
        self._pulses = pulses if pulses is not None else []
        self.baking_temp = baking_temp

        # Validate: if there are BAKE pulses, baking_temp must be set
        has_bake = any(p.pulse_type == "BAKE" for p in self._pulses)
        if has_bake and self.baking_temp is None:
            raise ValueError(
                "Scenario contains BAKE pulses but baking_temp was not set. "
                "Pass baking_temp=<value_in_K> to Scenario()."
            )

    @property
    def pulses(self) -> List[Pulse]:
        return self._pulses

    def to_txt_file(self, filename: str):
        df = pd.DataFrame(
            [
                {
                    "pulse_type": pulse.pulse_type,
                    "nb_pulses": pulse.nb_pulses,
                    "transition":pulse.transition,
                    "steady_state": pulse.steady_state,
                    "fraction_ss":pulse.fraction_ss,
                    "waiting": pulse.waiting,
                    "tritium_fraction": pulse.tritium_fraction,
                    "heat_scaling": pulse.heat_scaling,
                    "flux_scaling": pulse.flux_scaling,
                    "baking_temp": self.baking_temp,
                }
                for pulse in self.pulses
            ]
        )
        df.to_csv(filename, index=False)

    @staticmethod
    def from_txt_file(filename: str, old_format=False) -> "Scenario":
        if old_format:
            pulses = []
            with open(filename, "r") as f:
                for line in f:
                    # skip first line
                    if line.startswith("#"):
                        continue

                    # skip empty lines
                    if not line.strip():
                        continue

                    # assume this is the format
                    pulse_type, nb_pulses, ramp_up, steady_state, ramp_down, waiting = (
                        line.split()
                    )
                    pulses.append(
                        Pulse(
                            pulse_type=pulse_type,
                            nb_pulses=int(nb_pulses),
                            transition=[float(ramp_up),float(ramp_down)],
                            stead_STATE = [float(steady_state)],
                            steady_state=float(steady_state),
                            fraction_ss = [1],
                            waiting=float(waiting),
                        )
                    )
            return Scenario(pulses)
        df = pd.read_csv(filename)
        pulses = [
            Pulse(
                pulse_type=row["pulse_type"],
                nb_pulses=int(row["nb_pulses"]),
                transition=[float(row['ramp_up']),float(row['ramp_down'])],
                steady_state=float(row["steady_state"]),
                steady_STATE=[float(row["steady_state"])],
                fraction_ss=[1],
                waiting=float(row["waiting"]),
                tritium_fraction=float(row["tritium_fraction"]),
                heat_scaling=float(row.get("heat_scaling", 1.0)),
                flux_scaling=float(row.get("flux_scaling", 1.0)),
            )
            for _, row in df.iterrows()
        ]
        # Read baking_temp from CSV if present as a column (same value in all rows)
        if "baking_temp" in df.columns:
            baking_temp = float(df["baking_temp"].iloc[0])
        else:
            baking_temp = None  # will be validated in Scenario.__init__
        return Scenario(pulses, baking_temp=baking_temp)

    def get_row(self, t: float) -> int:
        """
        Returns the index of the pulse at time t.
        If t is greater than the maximum time in the scenario, a
        warning is raised and the last pulse index is returned.

        Args:
            t: the time in seconds

        Returns:
            the index of the pulse at time t
        """
        current_time = 0
        for i, pulse in enumerate(self.pulses):
            phase_duration = pulse.nb_pulses * pulse.total_duration
            if t < current_time + phase_duration:
                return i
            else:
                current_time += phase_duration

        warnings.warn(
            f"Time t {t} is out of bounds of the scenario file. Valid times are t < {self.get_maximum_time()}",
            UserWarning,
        )
        return i

    def get_pulse(self, t: float) -> Pulse:
        """
        Returns the pulse at time t.
        If t is greater than the maximum time in the scenario, a
        warning is raised and the last pulse is returned.

        Args:
            t: the time in seconds

        Returns:
            Pulse: the pulse at time t
        """
        row_idx = self.get_row(t)
        return self.pulses[row_idx]

    def get_pulse_type(self, t: float) -> str:
        """Returns the pulse type as a string at time t.

        Args:
            t: time in seconds

        Returns:
            pulse type (eg. FP, ICWC, RISP, GDC, BAKE)
        """
        return self.get_pulse(t).pulse_type

    def get_maximum_time(self) -> float:
        """Returns the maximum time of the scenario in seconds.

        Returns:
            the maximum time of the scenario in seconds
        """
        return sum([pulse.nb_pulses * pulse.total_duration for pulse in self.pulses])

    def get_time_start_current_pulse(self, t: float):
        """Returns the time (s) at which the current pulse started.

        Args:
            t: the time in seconds

        Returns:
            the time at which the current pulse started
        """
        pulse_index = self.get_row(t)
        return sum(
            [
                pulse.nb_pulses * pulse.total_duration
                for pulse in self.pulses[:pulse_index]
            ]
        )

    # TODO this is the same as get_time_start_current_pulse, remove
    def get_time_till_row(self, row: int) -> float:
        """Returns the time (s) until the row in the scenario file.

        Args:
            row: the row index in the scenario file

        Returns:
            the time until the row in the scenario file
        """
        return sum(
            [pulse.nb_pulses * pulse.total_duration for pulse in self.pulses[:row]]
        )

    # TODO remove
    def get_pulse_duration_no_waiting(self, row: int) -> float:
        """Returns the total duration (without the waiting time) of a pulse in seconds for a given row in the file.

        Args:
            row: the row index in the scenario file

        Returns:
            the total duration of the pulse in seconds
        """
        return self.pulses[row].duration_no_waiting

    # TODO remove
    def get_pulse_duration(self, row: int) -> float:
        """Returns the total duration of a pulse in seconds for a given row in the file.

        Args:
            row: the row index in the scenario file

        Returns:
            the total duration of the pulse in seconds
        """
        return self.pulses[row].total_duration
