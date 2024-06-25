import argparse
import pandas as pd
import numpy as np
from joblib import load
from sklearn.utils import resample

parser = argparse.ArgumentParser()
parser.add_argument("--m", type=str, help="Model Name")
parser.add_argument("--i", type=str, help="file name")

args = parser.parse_args()

modelName = "model"
if (args.m):
    modelName = str(args.m)

fileName = "data.csv"
if (args.i):
    fileName = str(args.i)

# Load the data
column_names = ['Ne', 'Emean_exhyt', 'Fix_index', 'Mlocus_homozegosity_mean', 'Mlocus_homozegosity_variance',
                    'Gametic_disequilibrium']
data = pd.read_csv(fileName, names = column_names, header = None, sep='\t')
X = np.array(data[['Emean_exhyt', 'Fix_index', 'Mlocus_homozegosity_mean', 'Mlocus_homozegosity_variance',
                      'Gametic_disequilibrium']])
y = np.array(data['Ne'])

# Load the model
model = load(modelName)

# bootstrap
n_iterations = 10000
n_size = int(len(data) * 0.50)
stats = list()
for i in range(n_iterations):
    X_sample, y_sample = resample(X, y, n_samples=n_size)
    res = model.predict(X_sample)
    stats.append(res)

# confidence intervals
stats = np.array(stats)  # 确保res是NumPy数组

# 计算95%置信区间
alpha = 0.975
p = ((1.0-alpha)/2.0) * 100
lower = np.percentile(stats, p)
p = (alpha+((1.0-alpha)/2.0)) * 100
upper = np.percentile(stats, p)
print('%.1f confidence interval %.1f and %.1f' % (alpha*100, lower, upper))