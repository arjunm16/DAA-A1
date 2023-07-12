#   Group Members
#   Arjun Muthiah       2019B3A70374H
#   Aryan Kapadia       2019B3A70412H
#   Harsh Vardhan Gupta 2019B3A70630H

import os
import subprocess
import pandas as pd

import math
import random
from typing import List, Tuple
from PIL import Image
from PIL import ImageDraw
import matplotlib.pyplot as plt

df = pd.DataFrame(columns = ["S.No.", "No. of Vertices", "No. of Notches", "Execution Time"])
# df.columns = ["S.No.", "No. of Vertices", "No. of Notches", "Execution Time"]

cpp_name = "htester_final_1.cpp"
num_ips = 5

if(os.system('g++ ' + cpp_name) == 0):
    for i in range(1, num_ips + 1):
        
        f = open("input" + str(i) + ".txt", "r")
        ln = f.readline().rstrip()
        f.close()

        proc = subprocess.run("./a.exe input"+ str(i) + ".txt", capture_output = True, text = True)
        
        notches, time = proc.stdout.rstrip().splitlines()
        notches = int(notches)
        time = float(time) / 1000
        
        #   [S.No., no. of verts, no. of notches, execution time]
        l = [i, int(ln), notches, time]
        df.loc[len(df)] = l

        # proc = subprocess.run("./a.exe ipL.txt", capture_output = True, text = True)
        # print(proc.stdout.rstrip())

        file1 = open("opinput" + str(i) + ".txt", "r")
        lines = file1.readlines()

        all_poly = []

        for line in lines:
            poly = []
            items = line.rstrip().split()
            for i in range(0, len(items),2):
                x = float(items[i])
                y = float(items[i+1])
                poly.append([x,y])
                
            # print(poly)
            all_poly.append(poly)
            
        # # print(all_poly)
        # print(all_poly)

        plt.figure()
        for i in range(len(all_poly)):
            all_poly[i].append(all_poly[i][0])
            xs, ys = zip(*(all_poly[i]))
            if i == (len(all_poly) - 1):
                plt.plot(xs,ys, color = "black")
            else:
                plt.plot(xs,ys, color = "red")
                # D:\Studies\4-2 Courses\CS F364 - DAA\Assignment\final version\figs
                # "D:/Studies/4-2 Courses/CS F364 - DAA/Assignment/final version/figs/img"
        # plt.savefig("/Users/aryankapadia/Desktop/2019B3A70630H_CSF364_A1/Ouputs/output15.jpg")
        plt.savefig("figs/img" + str(i) + ".jpg")
        # plt.show()

print(df)

# if(os.system('g++ ' + cpp_name) == 0):
#     ipfilename = "input1.txt"
#     f = open(ipfilename, "r")
#     ln = f.readline().rstrip()
#     f.close()

#     proc = subprocess.run("./a.exe "+ ipfilename, capture_output = True, text = True)
    
#     notches, time = proc.stdout.rstrip().splitlines()
#     notches = int(notches)
#     time = float(time) / 1000
    
#     #   [S.No., no. of verts, no. of notches, execution time]
#     l = [i, int(ln), notches, time]
#     df.loc[len(df)] = l

#     # proc = subprocess.run("./a.exe ipL.txt", capture_output = True, text = True)
#     # print(proc.stdout.rstrip())

#     file1 = open("op" + ipfilename, "r")
#     lines = file1.readlines()

#     all_poly = []

#     for line in lines:
#         poly = []
#         items = line.rstrip().split()
#         for i in range(0, len(items),2):
#             x = float(items[i])
#             y = float(items[i+1])
#             poly.append([x,y])
            
#         # print(poly)
#         all_poly.append(poly)
        
#     # # print(all_poly)
#     # print(all_poly)

#     plt.figure()
#     for i in range(len(all_poly)):
#         all_poly[i].append(all_poly[i][0])
#         xs, ys = zip(*(all_poly[i]))
#         if i == (len(all_poly) - 1):
#             plt.plot(xs,ys, color = "black")
#         else:
#             plt.plot(xs,ys, color = "red")
#             # D:\Studies\4-2 Courses\CS F364 - DAA\Assignment\final version\figs
#             # "D:/Studies/4-2 Courses/CS F364 - DAA/Assignment/final version/figs/img"
#     # plt.savefig("/Users/aryankapadia/Desktop/2019B3A70630H_CSF364_A1/Ouputs/output15.jpg")
#     plt.savefig("figs/img" + str(i) + ".jpg")
#     # plt.show()