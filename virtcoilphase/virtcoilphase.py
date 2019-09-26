'''Phase recon from multiple coil data using a virtual reference coil.
'''

import numpy as np
from scipy.signal import correlate2d
from tqdm import trange

def virtcoilphase(p, Cn=1, patch_size=(5, 5), coil_axis=-1):
    '''Absolute phase estimate via virtual reference coil.

    Parameters
    ----------
    p : array_like
        2D coil images in image space.
    Cn : array_like, optional
        Noise covariance matrix.
    patch_size : tuple, optional
        Size of patch to use to find most correlated point between
        all coils.
    coil_axis : int, optional
        Axis holding the coil data.  Other two axes will be assumed
        x and y.

    Returns
    -------
    phi_hat : array_like
        Absolute phase estimate.

    Notes
    -----
    Implements the algorithm described in [1]_.

    References
    ----------
    .. [1] Parker, Dennis L., et al. "Phase reconstruction from
           multiple coil data using a virtual reference coil."
           Magnetic resonance in medicine 72.2 (2014): 563-569.
    '''

    p = np.moveaxis(p, coil_axis, -1)
    sx, sy, nc = p.shape[:]

    # Thus absolute phase can be determined from multiple coil
    # measurements by:

    # (i) Creating a virtual reference coil that is sensitive over
    # the entire imaged volume using a complex weighted sum of the
    # individual coil measurements;

    # (a) the phase of the weights rotate the data from each coil
    # such that all coils have the same, hopefully zero, phase where
    # the coil sensitivities have maximum overlap;

    # Find a point where coil images are most correlated
    pabs = np.abs(p)
    pabs = pabs - np.mean(
        pabs.reshape((-1, nc)), axis=0)[None, None, :]
    corr = np.zeros((sx, sy, nc*nc))
    for ii in trange(nc, leave=False):
        for jj in range(nc):
            corr[..., jj*(ii+1)] = np.abs(correlate2d(
                pabs[..., ii], pabs[..., jj], mode='same'))
    corr = np.sqrt(np.sum(corr**2, axis=-1)) # is this okay?
    # corr = np.max(corr, axis=-1)
    x0 = np.unravel_index(np.argmax(corr.flatten()), corr.shape)
    # print(x0)

    # Now that we have the point, take the patch around it
    px, py = int(patch_size[0]/2), int(patch_size[1]/2)
    adjx, adjy = np.mod(patch_size, 2)
    pf = p[x0[0]-px:x0[0]+px+adjx, x0[1]-py:x0[1]+py+adjy, :]
    pf = np.mean(pf, axis=(0, 1)).reshape((-1, nc))
    phi_ref = -1*np.angle(
        np.sum(pf, axis=0)/np.sum(np.abs(pf), axis=0))
    # print(phi_ref)

    # (b) the weight magnitudes vary with the relative local SNR of
    # each coil;

    # w = np.zeros(p.shape, dtype=p.dtype)
    # den = np.sum(pabs, axis=-1)
    w = np.zeros(nc, dtype=p.dtype)
    den = np.sum(pabs.flatten())
    for jj in range(nc):
        # w[..., jj] = np.abs(p[..., jj])*np.exp(-1j*phi_ref[jj])/den
        w[jj] = np.abs(np.sum(pabs[..., jj].flatten()))*np.exp(
            -1*phi_ref[jj])/den

    # Now for the virtual coil
    # v = np.sum(w*p, axis=-1)
    v = p @ w

    # import matplotlib.pyplot as plt
    # plt.imshow(np.abs(v))
    # plt.show()

    # (ii) Replacing the individual coil phase with the virtual
    # reference coil phase while maintaining the object phase and
    # image noise of each coil by

    # (a) taking a low pass filter of the phase difference between
    # each coil and the reference coil;

    va = np.angle(v)
    win = np.hamming(sx)[:, None]*np.hamming(sy)[None, :]

    ax = (0, 1)
    Fv = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(
        np.angle(p) - va[..., None], axes=ax), axes=ax), axes=ax)
    dva = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(
        Fv*win[..., None], axes=ax), axes=ax), axes=ax)

    # (b) subtracting this low-pass-filtered phase difference
    # from the individual coil measurements;
    p = p*np.exp(-1j*dva)

    # (iii) using the inverse noise covariance to combine the phase
    # referenced complex measurements from all coils.

    # Two cases: noise covariance or no
    if isinstance(Cn, np.ndarray):
        iCn = np.linalg.pinv(Cn)
        val = np.zeros((sx, sy), dtype=p.dtype)
        for ii in range(nc):
            for jj in range(nc):
                val += p[..., ii]*iCn[ii, jj]*np.conj(p[..., jj])
        lamda = np.max(np.abs(val.flatten()))
        return va + np.angle(val/lamda)

    # No noise covariance is easy:
    val = np.sum(p*np.conj(p), axis=-1)
    lamda = np.max(np.abs(val.flatten()))
    return va + np.angle(val/lamda)
