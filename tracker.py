import pandas as pd
from datetime import datetime, timedelta
import argparse
import os

def read_simulation_log(file_path):
    with open(file_path, 'r') as file:
        reading = False
        data = []
        columns = None
        for line_number, line in enumerate(file, 1):
            if "Step Time Temp Press Enthalpy" in line and not reading:
                reading = True
                columns = line.strip().split()
                continue

            if reading:
                if "Step Time Temp Press Enthalpy" in line or "Loop time of" in line:
                    reading = False
                    continue
                else:
                    # 将每行数据拆分为一个列表，然后将其添加到数据中
                    data.append(list(map(float, line.strip().split())))

    # 创建一个DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df


def monitor_progress(start_time_str, target_steps, df):
    start_time = datetime.fromisoformat(start_time_str)
    current_time = datetime.now()
    elapsed_time = current_time - start_time

    # 计算当前进度
    current_step = df['Step'].iloc[-1]
    progress = min(current_step / target_steps * 100, 100)

    # 预估结束时间和所需剩余时间
    if progress < 100:
        estimated_end_time = start_time + timedelta(seconds=elapsed_time.total_seconds() * target_steps / current_step)
        estimated_remaining_time = estimated_end_time - current_time
    else:
        estimated_end_time = current_time
        estimated_remaining_time = timedelta(seconds=0)

    # 选择新DataFrame的行
    selected_rows = [0] + [int(i * len(df) / 4) for i in range(1, 4)] + [-1]
    if len(df) < 5:
        selected_rows = list(range(len(df)))
    new_df = df.iloc[selected_rows][['Step', 'Time', 'Temp', 'Press', 'Enthalpy', 'TotEng', 'Volume', 'Density']]

    print(f"Progress: {progress:.2f}%")
    print("*-" * 25)
    print(f"Start Time: {start_time}")
    print(f"Current Time: {current_time}")
    print(f"Estimated End Time: {estimated_end_time}")
    print(f"Estimated Remaining Time: {str(estimated_remaining_time).split('.')[0]}")
    print("*-" * 25)
    print("Progress Data:")
    col_widths = [10, 10, 10, 10, 10, 10, 10, 10]
    header = "".join([f"{col:<{col_widths[i]}}" for i, col in enumerate(new_df.columns)])
    print(header)
    for index, row in new_df.iterrows():
        row_str = "".join([f"{value:<{col_widths[i]}}" for i, value in enumerate(row)])
        print(row_str)


#extract basic simulation information from .in file
def extract_simulation_info(file_path):
    info = {
        "pair_style": None,
        "dielectric": None,
        "temperature": None,
        "timestep": None
    }

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        line = line.strip()

        # 提取力场样式
        if line.startswith('pair_style'):
            parts = line.split()
            info['pair_style'] = " ".join(parts[1:])

        # 提取介电常数
        elif line.startswith('dielectric'):
            parts = line.split()
            info['dielectric'] = float(parts[1])

        # 提取温度
        elif line.startswith('velocity') and 'create' in line:
            parts = line.split()
            info['temperature'] = float(parts[3])

        # 提取时间步长
        elif line.startswith('timestep'):
            parts = line.split()
            info['timestep'] = float(parts[1])

    return info

def print_info(info):
    print("Simulation Information:")
    print(f"Pair Style: {info['pair_style']} ")
    print(f"Dielectric Constant: {info['dielectric']} ")
    print(f"Temperature: {info['temperature']} K")
    print(f"Time Step: {info['timestep']} fs ")

#extract simulation step information from .in file
def extract_simulation_steps(file_path):
    steps = []
    total_steps = 0
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        line = line.strip()

        if line.startswith("fix") and ("nvt" in line or "npt" in line):
            parts = line.split()
            temperature_start = float(parts[5])
            temperature_end = float(parts[6])
            pressure = None
            if "iso" in parts:
                iso_index = parts.index("iso")
                pressure = (float(parts[iso_index + 1]), float(parts[iso_index + 2]))

            process = None
            if temperature_start != temperature_end:
                process = "Heating" if temperature_end > temperature_start else "Cooling"

            step = {
                "type": parts[3],
                "temperature_start": temperature_start,
                "temperature_end": temperature_end,
                "pressure": pressure,
                "process": process
            }
            steps.append(step)

        elif line.startswith("run"):
            parts = line.split()
            current_steps = int(parts[1])
            total_steps += current_steps
            steps[-1]['steps'] = current_steps

    return steps, total_steps


def print_steps(steps, total_steps):
    print("Simulation Steps:")
    for step in steps:
        print(f"Type: {step['type']}")
        print(f"Temperature: {step['temperature_start']} K to {step['temperature_end']} K")
        if step['pressure'] is not None:
            print(f"Pressure: {step['pressure'][0]} atm to {step['pressure'][1]} atm")
        if step['process']:
            print(f"Process: {step['process']}")
        print(f"Steps: {step['steps']} ({step['steps']} fs)")
        print("-" * 50)

    print(f"Total Steps: {total_steps} ({total_steps} fs)")



def main():
    parser = argparse.ArgumentParser(description="Simulation Monitor")
    parser.add_argument("log_file_path", help="Path to the log file")
    parser.add_argument("in_file_path", help="Path to the .in file")
    parser.add_argument("--info", action="store_true", help="Extract and print basic simulation information")
    parser.add_argument("--monitor", action="store_true", help="Monitor the progress")

    args = parser.parse_args()

    if args.info:
        info = extract_simulation_info(args.in_file_path)
        print_info(info)

    if args.monitor or args.info:
        steps, total_steps = extract_simulation_steps(args.in_file_path)
        if args.info:
            print_steps(steps, total_steps)

        if args.monitor:
            creation_time_stamp = os.path.getctime(args.in_file_path)
            creation_time = datetime.fromtimestamp(creation_time_stamp)
            start_time_str = creation_time.isoformat()

            df_in = read_simulation_log(args.log_file_path)
            monitor_progress(start_time_str, total_steps, df_in)

if __name__ == "__main__":
    main()


# python tracker.py /home/weilin/radonpy_file/radonpy_org/RadonPy-develop/sim_mix/work_dir_20230816_182311/radon_md.log /home/weilin/radonpy_file/radonpy_org/RadonPy-develop/sim_mix/work_dir_20230816_182311/radon_lmp.in --monitor
