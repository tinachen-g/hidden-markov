import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('posteriors.csv', header=None)

posterior_probs_h = data.iloc[0].values

sequence_positions = np.arange(1, len(posterior_probs_h) + 1)

plt.figure(figsize=(10, 6))
plt.plot(sequence_positions, posterior_probs_h, label='Posterior Probability for h', color='blue')

viterbi_intervals = [(66, 415), (527, 720), (951, 1000)]

for start, end in viterbi_intervals:
    plt.axvspan(start, end, color='red', alpha=0.3, label='Viterbi Interval' if start == viterbi_intervals[0][0] else "")

plt.xlabel('Sequence Position t')
plt.ylabel('P(π_t = h | x)')
plt.margins(x=0)
plt.ylim(0, 1.2)
plt.title('Posterior Probabilities for State h Over Sequence Position')
plt.legend()
plt.grid(True)

# Show plot
plt.show()


# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np


# data = pd.read_csv('posteriors.csv', header=None)  # Load the CSV without header
# posterior_probs = data.values.flatten()

# sequence_positions = np.arange(1, len(posterior_probs) + 1)

# plt.figure(figsize=(10, 6))
# plt.plot(sequence_positions, posterior_probs, label='Posterior Probabilities', color='blue')
# # 66,415
# # 527,720
# # 951,1000

# viterbi_intervals = [(66,415), (527,720),(951,1000) ] 

# for start, end in viterbi_intervals:
#     plt.axvspan(start, end, color='red', alpha=0.3, label='Viterbi Interval' if start == viterbi_intervals[0][0] else "")

# # Customize plot
# plt.xlabel('Sequence Position (t)')
# plt.ylabel('P(πt = h | x)')
# plt.title('Posterior Probabilities and Viterbi Intervals')
# plt.legend()
# plt.grid(True)

# # Show plot
# plt.show()


# import numpy as np
# import matplotlib.pyplot as plt

# # Step 1: Read the data
# # Assuming the CSV contains one line with probabilities and no headers
# posterior_probs = np.loadtxt('posteriors.csv', delimiter=',')

# # Create sequence positions (1-based indexing)
# sequence_positions = np.arange(1, len(posterior_probs) + 1)

# # Step 2: Plot posterior probabilities
# plt.figure(figsize=(10, 6))
# plt.plot(sequence_positions, posterior_probs, label='Posterior Probabilities', color='blue')

# viterbi_intervals = [(66,415), (527,720),(951,1000) ]  

# for start, end in viterbi_intervals:
#     plt.axvspan(start, end, color='red', alpha=0.3, label='Viterbi Interval' if start == viterbi_intervals[0][0] else "")

# # Customize plot
# plt.xlabel('Sequence Position (t)')
# plt.ylabel('P(πt = h | x)')
# plt.title('Posterior Probabilities and Viterbi Intervals')
# plt.legend()
# plt.grid(True)

# # Show plot
# plt.show()
