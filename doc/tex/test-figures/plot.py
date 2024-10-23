import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('exp-triple-dist.csv')
plt.hist(df)
plt.xlabel('Distance')
plt.ylabel('Frequency')
plt.savefig('exp-trip.pdf')
plt.show()

