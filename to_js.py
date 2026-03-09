import sys
import json

def convert_to_js(csv_file, js_file):
    steps = {}
    try:
        with open(csv_file, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) < 6: continue
                time = float(parts[0])
                # x, y, z are at index 3, 4, 5
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                
                time_str = f"{time:.3f}"
                if time_str not in steps:
                    steps[time_str] = []
                steps[time_str].append([x, y, z])
        
        # Write as a JS variable
        with open(js_file, 'w') as f:
            f.write("const simulationData = ")
            json.dump(steps, f)
            f.write(";")
        print(f"Converted {csv_file} to {js_file}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    convert_to_js('galaxy_step.csv', 'data.js')
