import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


k_propofol = 0.119 

k_dexmedetomidine = 0.00578 


PK_PARAMS = {
    'propofol': {
        'k': k_propofol,
        'halflife_min': np.log(2) / k_propofol,
        'vd_male': 75, 
        'vd_female': 65  
    },
    'dexmedetomidine': {
        'k': k_dexmedetomidine,
        'halflife_min': np.log(2) / k_dexmedetomidine,
        'vd_male': 91, 
        'vd_female': 91  
    }
}


def military_to_minutes(time_str):
    """
    Converts HHMM (e.g. '1302') to integer minutes after midnight.
    """
    if len(time_str) == 4:
        hours = int(time_str[:2])
        minutes = int(time_str[2:])
        return hours * 60 + minutes
    else:
        raise ValueError(f"Invalid time format: {time_str}")


def handle_midnight(start_m, end_m):
    """
    If end time < start time, assume it crossed midnight => end_m += 24*60.
    """
    if end_m < start_m:
        end_m += 24 * 60
    return end_m


def mg_to_mcg(mg):
    return mg * 1000.0


for med, params in PK_PARAMS.items():
    print(f"{med.capitalize()} half-life: {params['halflife_min']:.2f} minutes")


def is_steady_state(duration_minutes, halflife):
    """
    Return True if infusion ran for >= 4 half-lives.
    MODIFIED: Added halflife parameter to account for different medications
    """
    return duration_minutes >= 4 * halflife


def calc_bolus_concentration(t, t_bolus, dose_mcg, k, Vd):
    """
    Bolus at time t_bolus => (dose_mcg / Vd) * exp(-k*(t - t_bolus))
    t, t_bolus in minutes (absolute time)
    dose_mcg in micrograms, Vd in liters, k in 1/min.
    """
    if t < t_bolus:
        return 0.0
    return (dose_mcg / Vd) * np.exp(-k * (t - t_bolus))


def calc_infusion_concentration(t, t_start, t_end, rate, k, Vd):
    """
    Infusion from t_start to t_end with constant rate 'rate' (in mcg/kg/min, as in original).
    t, t_start, t_end in minutes (absolute times).
    Vd in liters, k in 1/min.
    """
    if t < t_start:
        return 0.0
    if t >= t_end:
        infusion_duration = t_end - t_start
        c_end = (rate / (k * Vd)) * (1 - np.exp(-k * infusion_duration))
        return c_end * np.exp(-k * (t - t_end))
    else:
        infusion_duration = t - t_start
        return (rate / (k * Vd)) * (1 - np.exp(-k * infusion_duration))


def calc_total_concentration(t, events, k, Vd):
    """
    Sum concentration contributions from all events at absolute time t (in minutes).
    Each event is either a bolus or infusion.
    """
    total = 0.0
    for e in events:
        if e['type'] == 'bolus':
            total += calc_bolus_concentration(t, e['time_minutes'], e['dose_mcg'], k, Vd)
        elif e['type'] == 'infusion':
            total += calc_infusion_concentration(t,
                                                 e['start_minutes'],
                                                 e['end_minutes'],
                                                 e['rate'],  # already in mcg/kg/min
                                                 k, Vd)
        else:
            raise ValueError(f"Unknown event type: {e['type']}")
    return total


def calc_total_concentration_by_medication(t, patient, med_name):
    """
    Calculate the total concentration for a specific medication at time t.
    """
    total = 0.0

    for med in patient['medications']:
        if med['name'] == med_name:
            k = med['k']
            Vd = med['Vd']

            for evt in med['events']:
                if evt['type'] == 'bolus':
                    total += calc_bolus_concentration(
                        t, evt['time_minutes'], evt['dose_mcg'], k, Vd)
                elif evt['type'] == 'infusion':
                    total += calc_infusion_concentration(
                        t, evt['start_minutes'], evt['end_minutes'],
                        evt['rate'], k, Vd)

            break

    return total


patients = [
    {
        'PatientID': 'ExamplePt',
        'Weight': 59.0,
        'Male': False,
        'medications': [
            {
                'name': 'propofol',
                'events': [
                    {
                        'type': 'infusion',
                        'start': '1120',
                        'end': '1136',
                        'rate_mcg_kg_min': 60
                    },
                    {
                        'type': 'infusion',
                        'start': '1136',
                        'end': '1157',
                        'rate_mcg_kg_min': 120
                    },
                    {
                        'type': 'infusion',
                        'start': '1157',
                        'end': '1232',
                        'rate_mcg_kg_min': 100
                    },
                    {
                        'type': 'infusion',
                        'start': '1232',
                        'end': '1242',
                        'rate_mcg_kg_min': 120
                    },
                    {
                        'type': 'infusion',
                        'start': '1242',
                        'end': '1302',
                        'rate_mcg_kg_min': 100
                    }
                ]
            }
        ]
    }
]



def preprocess_patient_events(p):
    """
    Convert a patient's events from military time to absolute minutes for all medications,
    handle midnight, and set appropriate PK parameters for each medication.
    Returns tuple (earliest_start, latest_stop) considering all medications.

    MODIFIED: Now handles multiple medication types with different parameters
    """
    earliest_start = None
    latest_stop = None

    for med in p['medications']:
        med_name = med['name']

        if p['Male']:
            med['Vd'] = PK_PARAMS[med_name]['vd_male']
        else:
            med['Vd'] = PK_PARAMS[med_name]['vd_female']

        med['k'] = PK_PARAMS[med_name]['k']

        for evt in med['events']:
            if evt['type'] == 'infusion':
                start_m = military_to_minutes(evt['start'])
                end_m = military_to_minutes(evt['end'])
                end_m = handle_midnight(start_m, end_m)
                evt['start_minutes'] = start_m
                evt['end_minutes'] = end_m

                if med_name == 'propofol':
                    evt['rate'] = evt['rate_mcg_kg_min']
                elif med_name == 'dexmedetomidine':
        
                    if 'rate_mcg_kg_min' not in evt:
                        evt['rate_mcg_kg_min'] = evt['rate_mcg_kg_hr'] / 60.0
                    evt['rate'] = evt['rate_mcg_kg_min']


                if earliest_start is None or start_m < earliest_start:
                    earliest_start = start_m
                if latest_stop is None or end_m > latest_stop:
                    latest_stop = end_m

            elif evt['type'] == 'bolus':
                t_m = military_to_minutes(evt['time'])
                evt['time_minutes'] = t_m
                evt['dose_mcg'] = mg_to_mcg(evt['dose_mg'])


    if earliest_start is None:
        earliest_start = 0
        latest_stop = 0

    return earliest_start, latest_stop


time_seconds = np.arange(0, 2000 + 10, 10)
data = {"Time (s)": time_seconds}


meta_data = {
    "Patient ID": [],
    "Medication": [],
    "Weight (kg)": [],
    "Gender": [],
    "Infusion Duration (min)": [],
    "Steady State Achieved": [],
    "Initial Rate": [], 
    "Final Rate": [],  
    "Calculated Initial Concentration (mcg/mL)": []
}


def calculate_infusion_duration_minutes(earliest, latest):
    return latest - earliest


for p in patients:
    pid = p['PatientID']
    weight = p['Weight']
    gender_str = "Male" if p['Male'] else "Female"

    earliest_start, latest_stop = preprocess_patient_events(p)
    duration_min = calculate_infusion_duration_minutes(earliest_start, latest_stop)

    for med_name in ['propofol', 'dexmedetomidine']:
        med_found = False
        for med in p['medications']:
            if med['name'] == med_name:
                med_found = True


                k = med['k']
                Vd = med['Vd']


                half_life = PK_PARAMS[med_name]['halflife_min']
                ss = "Yes" if is_steady_state(duration_min, half_life) else "No"

=
                infusion_events = [evt for evt in med['events'] if evt['type'] == 'infusion']
                if len(infusion_events) > 0:
                    if med_name == 'propofol':
                        R_initial = infusion_events[0]['rate_mcg_kg_min']
                        R_final = infusion_events[-1]['rate_mcg_kg_min']
                        initial_rate_str = f"{R_initial} mcg/kg/min"
                        final_rate_str = f"{R_final} mcg/kg/min"
                    else:  # dexmedetomidine
                        R_initial = infusion_events[0]['rate_mcg_kg_hr'] if 'rate_mcg_kg_hr' in infusion_events[0] else \
                        infusion_events[0]['rate_mcg_kg_min'] * 60
                        R_final = infusion_events[-1]['rate_mcg_kg_hr'] if 'rate_mcg_kg_hr' in infusion_events[-1] else \
                        infusion_events[-1]['rate_mcg_kg_min'] * 60
                        initial_rate_str = f"{R_initial} mcg/kg/hr"
                        final_rate_str = f"{R_final} mcg/kg/hr"
                else:
                    initial_rate_str = "0"
                    final_rate_str = "0"

                concentration_list = []
                for t_sec in time_seconds:
                    abs_time_min = latest_stop + (t_sec / 60.0)
                    c_val = calc_total_concentration_by_medication(abs_time_min, p, med_name)
                    concentration_list.append(c_val)

                column_name = f"{pid}_{med_name}"
                data[column_name] = concentration_list
                C0 = concentration_list[0]

                meta_data["Patient ID"].append(pid)
                meta_data["Medication"].append(med_name.capitalize())
                meta_data["Weight (kg)"].append(weight)
                meta_data["Gender"].append(gender_str)
                meta_data["Infusion Duration (min)"].append(duration_min)
                meta_data["Steady State Achieved"].append(ss)
                meta_data["Initial Rate"].append(initial_rate_str)
                meta_data["Final Rate"].append(final_rate_str)
                meta_data["Calculated Initial Concentration (mcg/mL)"].append(C0)
                break

        if not med_found:
            continue


time_df = pd.DataFrame({"Time (s)": time_seconds})


propofol_data = {"Time (s)": time_seconds}
dex_data = {"Time (s)": time_seconds}


propofol_meta = {k: [] for k in meta_data.keys()}
dex_meta = {k: [] for k in meta_data.keys()}


for i, med_type in enumerate(meta_data["Medication"]):
    pid = meta_data["Patient ID"][i]
    if med_type.lower() == "propofol":
        column_name = f"{pid}_{med_type.lower()}"
        if column_name in data:
            propofol_data[pid] = data[column_name]

        for k in meta_data.keys():
            propofol_meta[k].append(meta_data[k][i])
    elif med_type.lower() == "dexmedetomidine":
        column_name = f"{pid}_{med_type.lower()}"
        if column_name in data:
            dex_data[pid] = data[column_name]

        for k in meta_data.keys():
            dex_meta[k].append(meta_data[k][i])


propofol_df = pd.DataFrame(propofol_data)
dex_df = pd.DataFrame(dex_data)
propofol_meta_df = pd.DataFrame(propofol_meta)
dex_meta_df = pd.DataFrame(dex_meta)
                       
with pd.ExcelWriter("patient_sedation_decay.xlsx") as writer:
    propofol_df.to_excel(writer, sheet_name="Propofol Concentration", index=False)
    dex_df.to_excel(writer, sheet_name="Dexmedetomidine Concentration", index=False)
    propofol_meta_df.to_excel(writer, sheet_name="Propofol Metadata", index=False)
    dex_meta_df.to_excel(writer, sheet_name="Dexmedetomidine Metadata", index=False)

print("Excel file 'patient_sedation_decay.xlsx' has been saved with separate sheets for each medication.")


time_minutes = time_seconds / 60.0
