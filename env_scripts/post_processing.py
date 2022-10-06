import damask                             
import sys
import numpy as np 
import subprocess
import os


path = 'results/'

# Check whether the specified path exists or not
isExist = os.path.exists(path)
#if not create it 
if not isExist:
  
  os.makedirs(path)
  
#the result file .hdf5 is an argument of the function
result_file=str(sys.argv[1])
f = damask.Result(result_file)



simulated_stress=[]
simulated_strain=[]
stress_uniaxial=[]
strain_uniaxial=[]
#for all time increments, get the stress and strain 
for i in range(len(np.array(f.times))):
    print('increment', i)
    f_current = f.view('increments',f.increments[i])  
    sigma_2= (np.array(f_current.get('sigma'))) #get Cauchy stress
    #print(np.shape(sigma_2))
    if not np.shape(sigma_2):
        f_current.add_stress_Cauchy()
        f_current.add_strain()
        f_current.add_equivalent_Mises('sigma')
        f_current.add_equivalent_Mises('epsilon_V^0.0(F)')
    
    sigma_2= (np.array(f_current.get('sigma'))) #get Cauchy stress   
    epsilon_2= (np.array(f_current.get('epsilon_V^0.0(F)'))) #get von mises strain
    
    if  not np.shape(sigma_2):
    	print('time', f.times[i])
    	#print(sigma_2)
    elif not np.shape(epsilon_2):
    	print('time', f.times[i])
    	#print(epsilon_2)
    else:
        sigma_2=np.mean(sigma_2,0) #get average Cauchy stress on all grains
        simulated_stress.append(np.sqrt(1.0/2.0*((sigma_2[0,0]-sigma_2[1,1])**2+(sigma_2[1,1]-sigma_2[2,2])**2+(sigma_2[2,2]-sigma_2[0,0])**2)+3*(sigma_2[0,1]**2+sigma_2[1,2]**2+sigma_2[2,0]**2)))
     #get Von mises of the average Cauchy stress of all grains
        #get stress in the xx direction 
        stress_uniaxial.append(sigma_2[0,0])
        epsilon_2= (np.array(f_current.get('epsilon_V^0.0(F)'))) #get strai
        
        epsilon_2=np.mean(epsilon_2,0) #get average strains on all grains
        strain_uniaxial.append(epsilon_2[0,0])
        simulated_strain.append(np.sqrt(1.0/2.0*((sigma_2[0,0]-sigma_2[1,1])**2+(sigma_2[1,1]-sigma_2[2,2])**2+(sigma_2[2,2]-sigma_2[0,0])**2)+3*(sigma_2[0,1]**2+sigma_2[1,2]**2+sigma_2[2,0]**2)))


dirname = os.getcwd()
#write stress strain data in result directory
filename = "results/stress_strain"
pathname = os.path.abspath(os.path.join(dirname, filename))
if pathname.startswith(dirname):
   if os.path.exists(pathname):
       os.remove(pathname)
   
   
#write stress strain data in result directory
np.savetxt('results/stress_strain', np.c_[ strain_uniaxial, stress_uniaxial ])
#np.savetxt('stress_strain.txt' , [strain_uniaxial, stress_uniaxial], fmt=['%f','%f'])


for fname in os.listdir(dirname):
    if fname.startswith("results/final_"):
        os.remove(os.path.join(dirname, fname))
    if fname.startswith("results/xi_"):
        os.remove(os.path.join(dirname, fname))
    if fname.startswith("results/c_vacancy"):
        os.remove(os.path.join(dirname, fname))
    if fname.startswith("results/precipitate_distribution"):
        os.remove(os.path.join(dirname, fname))




for k in range(len(np.array(f.times))):

	#get some variables and save them in textfiles
	f_last = f.view('increments',f.increments[k])
	final_radius=f_last.get('avg_precipitate_radius')
	final_von_mises_strain=f_last.get('epsilon_V^0.0(F)_vM')
	final_von_mises_stress=f_last.get('sigma_vM')
	final_vf=f_last.get('precipitate_volume_frac')
	final_xi=f_last.get('xi_sl')
	final_vacancy_concentration = f_last.get('c_vacancy')
	#each of the textfile corresponds to one time increment and prints the variable in all voxels of the microstructure
	np.savetxt('results/vf_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_vf])
	np.savetxt('results/vm_strain_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_von_mises_strain])
	np.savetxt('results/vm_stress_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_von_mises_stress])
	np.savetxt('results/radius_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_radius])
	np.savetxt('results/xi_sl_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_xi])
	np.savetxt('results/c_vacancy_{t}'.format(t=int(np.round(f.times[k],0))), np.c_[final_vacancy_concentration])




f_first = f.view('increments',f.increments[0]) 
#print(f_first.get('precipitate_density').shape)
precipitate_distribution = np.mean(f_first.get('precipitate_density'),0)
np.savetxt('results/precipitate_distribution_inc_{i}.txt'.format(i=0) , [precipitate_distribution])





f_last = f.view('increments',f.increments[-1]) 
precipitate_distribution = np.mean(f_last.get('precipitate_density'),0)
np.savetxt('results/precipitate_distribution_last.txt'.format(i=0) , [precipitate_distribution])
#np.savetxt('stress_strain.txt' , [strain_uniaxial, stress_uniaxial], fmt=['%f','%f'])

#save information in vtr files that can be displayed with paraview
for i in range(len(np.array(f.times))):
	f_current= f.view('increments',f.increments[i])
	f_current.save_VTK(output=[ 'avg_precipitate_radius', 'c_vacancy' , 'precipitate_volume_frac', 'total_precipitate_density', 'sigma_vM','epsilon_V^0.0(F)_vM', 'xi_sl'  ])   # Write defined datasets into .vtr for each increment
