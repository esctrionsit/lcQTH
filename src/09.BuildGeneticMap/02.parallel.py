from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED
import concurrent.futures
from time import sleep
import traceback
import sys
import os

MAXTHR = int(sys.argv[1])

with open("chrlen.txt") as f:
    lines = f.readlines()

def task(chridx):
    try:
        shell = "Rscript build.parallel.R GM/ " + str(chridx+1)
        os.system(shell)
    except Exception as e:
        #直接打印错误
        traceback.print_exc()

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=MAXTHR) as t: 
        all_task = []
        for i in range(len(lines)):
            all_task.append(t.submit(task, i))
        wait(all_task, return_when=ALL_COMPLETED)