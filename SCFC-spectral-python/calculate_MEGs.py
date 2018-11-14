import time
import spectral_graph as sg
import util

hcp_dir = util.get_absolute_path('./data/')

freqresp, ev, Vv, freqresp_out, FCmodel = sg.network_transfer_function(C, Ddk_conn, w=0)

evec = []
Vvec =[]
fqall = []
freq_model = []

start = time.time()
for i in fvec:
    w = 2*np.pi*i
    fq, ev, Vv, freqresp_out, _ = NetworkTransferFunction(C, Ddk_conn, w)
    fqall.append(fq)
    evec.append(ev)
    Vvec.append(Vv)
    freq_model.append(freqresp_out)

fqall=np.asarray(fqall)
evec = np.asarray(evec)
Vvec = np.asarray(Vvec)
freq_model = np.asarray(freq_model)
stop = time.time()
duration = stop - start
print('Computation time = ', duration)
np.moveaxis(Vvec,1,0).shape #86x40x57 just like matlab

# Plotting eigen values of each eigen vector
#mpl.figure(num=5)
fig_ev, ax_ev = mpl.subplots(1,3, figsize = (20,3), sharey = True)
for i in np.arange(0,evec.shape[1]):
    ax_ev[0].scatter(np.abs(evec[:,i]), np.ones(evec.shape[0])*(i+1), c = np.arange(0,fvec.size))

ax_ev[0].grid(True)
ax_ev[0].set_xlabel('Eigen Value')
ax_ev[0].set_ylabel('Eigen Vector #')
ax_ev[0].set_title('Corresponding Eigen Values')
#mpl.savefig('ev_Vv_plot.png', dpi = 300, format = 'png')

# Plotting eigen vector's frequency response with default parameters
ev_freqresp = np.abs(np.transpose(fqall))
#fig_ev, ax_ev = mpl.subplots()
ax_ev[1].imshow(mag2db(np.abs(ev_freqresp)), extent = [0, 40, 0, 57], aspect='auto')
#fig_ev.colorbar(ax1)
ax_ev[1].set_title('Freq. Response of Eigen Vectors - Default Parameters')
ax_ev[1].set_xlabel('Frequency (Hz)')
#mpl.savefig('freqresp_eigs_default.png', dpi = 300, format = 'png')

# eig vs. freq figure with set of optimized parameters
evec = []
Vvec =[]
fqall = []
freq_model = []

# opparam = loadmat('/home/axiezai/lab/SCFC_eeg/data/SCFC_opparam_HCP.mat')
# optau_e = opparam['output']['param'][0,1][0]
# optau_i = opparam['output']['param'][0,1][1]
# opalpha = opparam['output']['param'][0,1][2]
# opspeed = opparam['output']['param'][0,1][3]
# opgei = opparam['output']['param'][0,1][4]
# opgii = opparam['output']['param'][0,1][5]
# optauC = opparam['output']['param'][0,1][6]

for i in fvec:
    w = 2*np.pi*i
    fq, ev, Vv, freqresp_out, _ = NetworkTransferFunction(C, Ddk_conn, w) #,tau_e = optau_e, tau_i = optau_i, alpha = opalpha,
                                                          #speed = opspeed, gei = opgei, gii = opgii, tauC = optauC)
    fqall.append(fq)
    evec.append(ev)
    Vvec.append(Vv)
    freq_model.append(freqresp_out)

fqall=np.asarray(fqall)
evec = np.asarray(evec)
Vvec = np.asarray(Vvec)
freq_model = np.asarray(freq_model)
freq_model = np.transpose(freq_model)
np.moveaxis(Vvec,1,0).shape #86x40x57 just like matlab

# Plotting eigen vector's frequency response
ev_freqresp = np.abs(np.transpose(fqall))
ax1 = ax_ev[2].imshow(mag2db(np.abs(ev_freqresp)), extent = [0, 40, 0, 57], aspect='auto')
fig_ev.colorbar(ax1)
ax_ev[2].set_title('Freq. Response of Eigen Vectors - Optimized Parameters')
ax_ev[2].set_xlabel('Frequency (Hz)')
mpl.savefig('freqresp_eigs.png', dpi = 300, format = 'png')
