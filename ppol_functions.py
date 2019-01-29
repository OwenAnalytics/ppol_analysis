def plot_waveforms(trz, trn, tre, tp, nb_event, taille, my_cmap, path_to_fig):

    dt = np.abs(tp - trz.stats.starttime)  # Time between stream starrtime and arrivaltime

    # Plot filtered waveform
    f, axarr = plt.subplots(3, sharex = true)
    axarr[0].plot(trz.times(), trz.data, color = my_cmap.colors[250])
    axarr[0].set_ylabel('%s' %trz.stats.channel)
    axarr[0].axvline(dt, ymin = min(trz.data), ymax = max(trz.data), color = my_cmap.colors[20])
    axarr[0].set_title('filtered waveform \n %s - %s' %(trz.stats.starttime, trz.stats.endtime))
    axarr[1].plot(trn.times(), trn.data, color = my_cmap.colors[250])
    axarr[1].set_ylabel('%s' %trn.stats.channel)
    axarr[1].axvline(dt, ymin = min(trz.data), ymax = max(trz.data), color = my_cmap.colors[20])
    axarr[2].plot(tre.times(), tre.data, color = my_cmap.colors[250])
    axarr[2].set_ylabel('%s' %tre.stats.channel)
    axarr[2].axvline(dt, ymin = min(trz.data), ymax = max(trz.data), color = my_cmap.colors[20])
    axarr[2].set_xlabel('time (s)')
    for i in range(len(axarr)):
            axarr[i].set_xlim(0, trz.stats.npts * trz.stats.delta )

    f.savefig(path_to_fig+'/lsb%s_figs/event%03d_filt' %(ns, nb_event), format = 'ps')
    plt.show()

    # Plot chosen analysis window
    f, axarr = plt.subplots(3, sharex = true)

    for i in range(len(axarr)):
            axarr[i].set_xlim(dt - taille, dt + 2*taille)
            lims = axarr[i].get_xlim()
            we = np.where( (trz.times()  > lims[0]) & (trz.times() < lims[1]))
            lim = np.max(np.abs(trz.data[we]))
            axarr[i].set_ylim( -lim, lim)
            #axarr[i].set_xlim(dt - 2*taille, dt + 5*taille)
            #axarr[i].set_ylim( tr_z.data[we].min(), tr_z.data[i].max())

    axarr[0].plot(trz.times(), trz.data, color = my_cmap.colors[250])
    axarr[0].set_ylabel('%s' %trz.stats.channel)
    axarr[0].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[0].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
    axarr[0].set_title('analysis window of filtered waveform')

    axarr[1].plot(trn.times(), trn.data, color = my_cmap.colors[250])
    axarr[1].set_ylabel('%s' %tr_n.stats.channel)
    axarr[1].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[1].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])

    axarr[2].plot(tre.times(), tre.data, color = my_cmap.colors[250])
    axarr[2].set_ylabel('%s' %tre.stats.channel)
    axarr[2].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[2].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
    axarr[2].set_xlabel('time (s)')

    f.savefig(path_to_fig +'/lsb%s_figs/event%03d_window' %(ns, nb_event), format = 'ps')
    plt.show()

    # Adjust analysis window if necessary
    #ch = 0
    ch = float(input('change in arrival time?'))

    f, axarr = plt.subplots(3, sharex = true)

    for i in range(len(axarr)):
            axarr[i].set_xlim(dt - taille, dt + 2*taille)
            lims = axarr[i].get_xlim()
            we = np.where( (trz.times() > lims[0]) & (trz.times() < lims[1]))
            lim = np.max(np.abs(trz.data[we]))
            axarr[i].set_ylim( -lim, lim)
            #axarr[i].set_xlim(dt - 2*taille, dt + 5*taille)
            #axarr[i].set_ylim( tr_z.data[we].min(), tr_z.data[i].max())

    axarr[0].plot(trz.times(), trz.data, color = my_cmap.colors[250])
    axarr[0].set_ylabel('%s' %tr_z.stats.channel)
    axarr[0].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[0].axvline(dt + ch + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
    axarr[0].set_title('analysis window of filtered waveform')

    axarr[1].plot(trn.times(), trn.data, color = my_cmap.colors[250])
    axarr[1].set_ylabel('%s' %tr_n.stats.channel)
    axarr[1].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[1].axvline(dt + ch +taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])

    axarr[2].plot(tre.times(), tre.data, color = my_cmap.colors[250])
    axarr[2].set_ylabel('%s' %tr_e.stats.channel)
    axarr[2].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
    axarr[2].axvline(dt + ch + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
    axarr[2].set_xlabel('time (s)')

    
    f.savefig(path_to_fig+'/lsb%s_figs/event%03d_window' %(ns, nb_event), format = 'ps')
    plt.show()
    return ch

   
def ppol_analysis(tre, trn, trz, my_cmap)
    # Put traces in mean centered e,n,z matrix
    mat = np.column_stack((tre.data, trn.data, trz.data)) 

    for j in range(3):
            mat[:,j] = mat[:,j] - np.mean(mat[:,j])

    # Calculate covariance
    cov = np.dot((np.transpose(mat)), mat) * (1/trim_z.stats.npts)

    # Eigenvalues and vectors
    val, vec = np.linalg.eig(cov)


    args = np.argsort(val)[::-1]  # Get sorted args
    val = np.sort(val)[::-1]  # Sort eigenvalues
    vecs = np.column_stack((vec[:,args[0]],vec[:,args[1]],\
            vec[:,args[2]]))  # Sort eigenvectors according to eigenvalues

    # Normalize first eigenvector - should already be done by numpy
    vec1 = vecs[:,0]
    #nvec1 =  np.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2) 
    #vec1 = vec1 / nvec1

    # Calculate angle of polarization
    if vec1[0] <= 0 and vec1[1] < 0:
        Theta = (np.pi + np.arctan(vec1[0] / vec1[1])) * 180 / np.pi
    elif vec1[0] < 0  and vec1[1] > 0:
        Theta = ( 2 * np.pi + np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi
    elif vec1[0] >= 0 and vec1[1] > 0:
        Theta = (np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi
    elif vec1[0] >= 0 and vec1[1] <  0:
        Theta = ( np.pi + np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi

    # Degree of rectilinarity 3D

    rect = 1 - (val[1] + val[2]) / (2*val[0])

    # Horizontal covariance matrix
    cov_h = np.zeros((2,2))
    cov_h[0,0] = cov[0,0]
    cov_h[0,1] = cov[0,1]
    cov_h[1,0] = cov[1,0]
    cov_h[1,1] = cov[1,1]

    valh, vech = np.linalg.eig(cov_h)  # Eigenvalues and eigenvectors

    # Degree of horizontal rectilinarity
    CpH = 1 - np.min(valh)/np.max(valh)

    # 2D azimuth
    argsh = np.argsort(valh)[::-1]  # Get sorted args
    valh = np.sort(valh)[::-1] # Sort eigenvalues
    vecsh = np.column_stack((vech[:,argsh[0]],vech[:,argsh[1]]))  # Sort eigenvectors according to eigenvalues

    # Normalize first eigenvector - should already be done by numpy

    vec1h = vecsh[:,0]
    
    if vec1h[0] >= 0 and vec1h[1] > 0:
        ThetaH = (np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
    elif vec1h[0] >= 0  and vec1h[1] < 0:
        ThetaH = (np.pi + np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
    elif vec1h[0] < 0 and vec1h[1] > 0:
        ThetaH = (2 * np.pi + np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
    elif vec1h[0] < 0 and vec1h[1] <  0:
        ThetaH = ( np.pi + np.arctan(vec1[0] / vec1[1])) * 180 / np.pi
    

    # Resolution of the ambiguity between thetaH and theta
    ep1 = Theta - ThetaH
    ep2 = ThetaH - Theta

    if ep1 > 140:
        Theta = Theta - 180

    if ep2 > 140:
        Theta = Theta + 180

    # Signal to noise ration
    SNR = (np.max(valh) - np.min(valh)) / np.min(valh)

    # Angle of vertical polarization

    # Longitudinatl component
    dnorthr = (ThetaH + 180.) * np.pi / 180
    L = np.cos(dnorthr) * (trim_n.data) + np.sin(dnorthr) * (trim_e.data)

    f, axarr = plt.subplots(2)
    # Horizontal particle motion
    axarr[0].plot(tre.data - np.mean(tre.data), trn.data - np.mean(trn.data), color = my_cmap.colors[250])
    axarr[0].set_title('Horizontal particle motion')
    axarr[0].set_ylabel('BHN')
    axarr[0].set_xlabel('BHE')
    limsx = axarr[0].get_xlim()
    limsy = axarr[0].get_ylim()
    limxy = np.max([np.max(np.abs(limsx)), np.max(np.abs(limsy))])
    axarr[0].set_ylim(-limxy, limxy)
    axarr[0].set_xlim(-limxy, limxy)
    x_lin = np.linspace(-limxy, limxy, 10)
    m = vec1h[1] / vec1h[0]
    axarr[0].plot(x_lin, m * x_lin, color = my_cmap.colors[50])
    

    # Vertical particle motion
    axarr[1].plot(L, trz.data, color = my_cmap.colors[250])
    axarr[1].set_title('Vertical particle motion')
    axarr[1].set_ylabel('BHZ')
    axarr[1].set_xlabel('L')
    
    limsx = axarr[1].get_xlim()
    limsy = axarr[1].get_ylim()
    limxy = np.max([np.max(np.abs(limsx)), np.max(np.abs(limsy))])
    axarr[1].set_ylim(-limxy, limxy)
    axarr[1].set_xlim(-limxy, limxy)
    # Set subfigures in quadratic form
    for ax in axarr:
        ax.set(adjustable='box-forced', aspect='equal')
    plt.tight_layout()
    plt.show()
    f.savefig(path_to_fig + '/Event%03d_motion' %(diri, ns, ne), format = 'PS')

    Vpol = np.arccos(vec1[2]) * 180 / np.pi

    if Vpol <= 90 and Vpol >= 0:
        Vpol = Vpol
    elif vec1[0] < 0 and vec1[1] < 0  and vec1[2] > 0 and Vpol < 90:
        Vpol = 180 - Vpol;
    elif vec1[0] > 0 and vec1[1] > 0 and vec1[2] > 0 and Vpol < 90:
        Vpol = 180 - Vpol
    elif vec1[0] < 0 and vec1[1] > 0 and vec1[2] > 0 and Vpol < 90:
        Vpol = 180 - Vpol


    if Vpol <= 90 and Vpol >= 0:
        Vpol = -1 * ( Vpol - 180 )
        Theta = Theta - 180
        ThetaH = ThetaH - 180

    VV = np.column_stack((L - np.mean(L), trim_z.data - np.mean(trim_z.data)))
    cov_v = np.dot(np.transpose(VV), VV) * (1/trim_z.stats.npts)

    valv, vecv = np.linalg.eig(cov_v)

    CpZ = 1 - np.min(valv) / np.max(valv)

    err2d = np.rad2deg(np.arctan(np.sqrt(np.min(valh) / np.max(valh)))) #Calculate error
    errV = np.rad2deg(np.arctan(np.sqrt(np.min(valv) / np.max(valv))))
    err3d = np.rad2deg(np.arctan(np.sqrt(val[2] / (val[1] + val[0]))))

    if Theta < 0:
        Theta += 360
    if ThetaH < 0:
        ThetaH += 360

    #print('BAZmeas3D %s' %Theta)
    #print('BAZmeas2D %s' %ThetaH)
    #print('CpH %s' %CpH)
    #print('CpZ %s' %CpZ)
    #print('ErH2D %s' %err2d)
    #print('ErH3D %s' %err3d)
    plt.close('all')
    
    return SNR, Theta, ThetaH, CpH, CpZ, err3d, err2d, baz


    

