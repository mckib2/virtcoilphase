'''Simple example showing proof of concept.'''

# import numpy as np
import matplotlib.pyplot as plt
from bart import bart # pylint: disable=E0401

from phaserecon import phaserecon

if __name__ == '__main__':

    imspace = bart(1, 'phantom -x64 -s8').squeeze()
    phi_hat = phaserecon(imspace)

    plt.imshow(phi_hat)
    plt.show()
