import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    print("CL: Start with time plotting.")
    time_report = sys.argv[1]
    if (time_report == None):
        print("Finishing")
        return

    time_input = pd.read_csv(time_report, sep= '\t\t')
    
    plt.figure(figsize=(8,4))
    plt.bar(time_input["Sample name"], time_input["Total time"])

    plt.xlabel("Sample name")
    plt.ylabel("Execution time (s)")
    plt.title("Pipeline execution time for each sample")

    plt.tight_layout()
    plt.show()

    plt.savefig("./Output/INT26_TIME_REPORT_GRAPH.png", dpi=300)

    print("CL: Done with time plotting.")

if __name__ == "__main__" : main()