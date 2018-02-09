import numpy as np
import matplotlib.pyplot as plt
import mode_analysis_code
import coldatoms
import time
import datetime
# Ek = np.load('kinetic_energy.npy')
# Ep_c = np.load('potential_energy_coulomb.npy')
# Ep_t = np.load('potential_energy_trap.npy')
# Ep = np.load('potential_energy.npy')
Etotal = np.load('total_energy.npy')

# Etotal_part = Etotal[0:Etotal.shape[0]//500]
# fig=plt.figure()
# plt.plot(Etotal_part[0:Etotal_part.shape[0]:1]*1e19, 'r', ms=0.05)
# plt.ylim([-1.7782, -1.7786])
# plt.show()
#plt.savefig('energy_evolution.pdf')

# trajectories=np.load('trajectories.npy')
# velocities=np.load('velocities.npy')
#
# epslion0 = 8.854187e-12
# q = 1.602176565e-19
# #import from mode_analysis_code
# kz=2.0 * mode_analysis.Coeff[2]
# delta=mode_analysis.Cw
# omega=mode_analysis.wrot
# phi0=np.pi/2.0
# #print(trajectories.shape[0],trajectories.shape[1],trajectories.shape[2])
# def coulomb_potential(q,epslion,trajectories):
#     potential=np.zeros(trajectories.shape[0])
#     N = trajectories.shape[1]
#     K=q*q/(4*np.pi*epslion)
#
#     for step in range(trajectories.shape[0]):
#         position_difference=np.zeros((int(N*(N-1)/2),3))
#         distance=np.zeros(int(N*(N-1)/2))
#         j = 0
#         for n in range(N):
#             for n_behind in range(N-n-1):
#                 position_difference[j,:]=trajectories[step,n,:]-trajectories[step,n+n_behind+1,:]
#                 j+=1
#         distance = np.sqrt(np.sum(np.square(position_difference),axis=1))
#         potential[step] =np.sum(np.divide(K,distance))
#     return potential
# def penning_trap_potential(q,kz,delta,omega,phi0,trajectories):
#     N = trajectories.shape[1]
#     potential_series=np.zeros((trajectories.shape[0],trajectories.shape[1]))
#     potential=np.zeros(trajectories.shape[0])
#     #define potential parameters
#     kx = (0.5+delta)*kz
#     ky = (0.5-delta)*kz
#     phi= phi0
#     dt = 1.0e-9
#     sampling_period = 2.5e-7
#     delta_t = dt*sampling_period
#
#     for step in range(trajectories.shape[0]):
#         x=trajectories[step,:,0]
#         y=trajectories[step,:,1]
#         z=trajectories[step,:,2]
#         s=np.sin(phi)
#         c=np.cos(phi)
#         potential_series[step]=0.5*kz*np.square(z)-0.5*(kx*np.square(c*x-s*y)+ky*np.square(s*x+c*y))
#         phi+=delta_t*omega
#     potential_series*=q
#     potential=np.sum(potential_series,axis=1)
#     return potential
#
# kinetic_energy = 0.5*m_Be*np.sum(np.sum(np.square(velocities),axis=2),axis=1)
# potential_energy_coulomb = coulomb_potential(q,epslion0,trajectories)
# potential_energy_trap =penning_trap_potential(q,kz, delta, omega, phi0,trajectories)
# potential_energy = potential_energy_coulomb+potential_energy_trap
# total_energy = kinetic_energy+potential_energy
#
# np.save('kinetic_energy.npy',kinetic_energy)
# np.save('potential_energy_coulomb.npy',potential_energy_coulomb)
# np.save('potential_energy_trap.npy',potential_energy_trap)
# np.save('potential_energy.npy',potential_energy)
# np.save('total_energy.npy',total_energy)
#
# def dev(array):
#     return np.sqrt(np.var(array))/np.average(array)


## calculation energy 2spectrum
print(datetime.datetime.now())
psd = np.abs((np.fft.fft(Etotal)/Etotal.shape[0])**2)
psd = psd[0:psd.size//2:1] + psd[2*(psd.size//2):psd.size//2:-1]
print(datetime.datetime.now())

sampling_period= 2.5e-7
nu_nyquist = 0.5 / sampling_period
nu_axis = np.linspace(0.0, nu_nyquist, Etotal.shape[0] // 2)

fig = plt.figure()
plt.semilogy(nu_axis / 1.0e6, psd,
             linewidth=0.75, color='blue', zorder=-1)
plt.xlabel(r'$\nu / \rm{MHz}$')
plt.ylabel(r'PSD($z$)')
plt.show()
# plt.savefig('energy_spectrum.pdf')
plt.xlim([0.0, 0.1])
plt.ylim([1.0e-17, 1.0e-13])
