from mc import *
from tqdm import tqdm

N = int(1e6)

def generate_random_3d_vector(length):
    # Generate a random 3D unit vector
    theta = 2 * np.pi * np.random.random()  # azimuthal angle
    phi = np.arccos(2 * np.random.random() - 1)  # polar angle
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    # Convert the unit vector to a vector of the specified length
    vector = length * np.array([x, y, z])
    return vector

# for _ in tqdm(range(int(N))):
#     generate_random_3d_vector(0.3)
    
for _ in tqdm(range(int(2))):
    theta = 2 * np.pi * np.random.random(N)  # azimuthal angle
    phi = np.arccos(2 * np.random.random(N) - 1)  # polar angle
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)

    # Convert the unit vector to a vector of the specified length
    vector = 0.3 * np.array([x, y, z]).T
    print(vector[0])