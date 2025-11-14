import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from itertools import combinations
import time

start_time = time.time()

def calculate_combination_means(data):
    num_seeds = len(data)-1
    combs = combinations(data, num_seeds)
    means = [sum(comb) / num_seeds for comb in combs]    
    return means

def corr_func(data, corr_len):
    print("data:", data.shape)
    
    data_1 = data[:, 0]  #xaxis
    data_2 = data[:, 1]  #yaxis
    data_3 = data[:, 2]  #zaxis

    print("data_1:", data_1.shape)

    result_x = np.correlate(data_1, data_1, mode='full')
    result_y = np.correlate(data_2, data_2, mode='full')
    result_z = np.correlate(data_3, data_3, mode='full')

    print("result_x:", result_x.shape)

    r1_x = result_x[result_x.size//2:]  
    r1_y = result_y[result_y.size//2:]  
    r1_z = result_z[result_z.size//2:]    

    print("r1_x:", r1_x.shape)

    corr_len_divide = np.arange(r1_x.shape[0], 0, -1)
    
    print("corr_len_divide:", corr_len_divide.shape)
    
    r1_x_norm = r1_x/corr_len_divide
    r1_y_norm = r1_y/corr_len_divide
    r1_z_norm = r1_z/corr_len_divide

    r2_x = r1_x_norm[:corr_len+1]
    r2_y = r1_y_norm[:corr_len+1]
    r2_z = r1_z_norm[:corr_len+1]
    
    return r2_x, r2_y, r2_z

temp_range = [300]
sub_path = ["./GK/0.1AlGaN"] 

all_seed_ranges = {
    "./GK/0.1AlGaN": np.asarray([2314, 3182, 3534, 4595, 4694, 6729]),
    } 


vol_array = np.asarray([91713.687]) # for alloys

sample_freq = 5

ev2J = 1.60218e-19
kB = 1.38e-23 #(J/K)
p2s = 1e-12
A2m = 1e-10
deltaT = 0.001*p2s

for ind_dir, sdir in enumerate(sub_path):
    V = vol_array[ind_dir] * A2m**3   
    for temp_ind, temp in enumerate(temp_range):  
        seed_range = all_seed_ranges[f"{sdir}"]   
         
        T = temp #K

        Flux_conv = (ev2J)**2 * (A2m/p2s)**2  
        Unit_conv = Flux_conv/kB/T/T/V*deltaT*sample_freq
        
        num_seeds = len(seed_range)
        corr_len = 2000 - 1
        corrfull_x = np.zeros([corr_len+1, num_seeds])
        corrfull_y = np.zeros([corr_len+1, num_seeds])
        corrfull_z = np.zeros([corr_len+1, num_seeds])
        file_path = f"./{sdir}/"
        print("Currently running:", temp, "K in path:", file_path, "with volume of:", vol_array[ind_dir])

        for seed_ind, seed in enumerate(seed_range):
            print("seed:", seed)
            raw_flux = np.loadtxt(file_path + "raw_heatcurrent_" + str(seed) + ".txt")
            raw_flux = raw_flux[:, 1:] * vol_array[temp_ind] # to remove time step column 
            corrfull_x[:, seed_ind], corrfull_y[:, seed_ind], corrfull_z[:, seed_ind] = corr_func(raw_flux, corr_len)    
        
        print("Temp:", temp, "end of check 1:", (time.time() - start_time))

        print("Temp:", temp, "data after corr in x, y and z:", corrfull_x.shape, corrfull_y.shape, corrfull_z.shape)
        
        np.save(file_path + "corr_before_mean_X_" + str(T) + "_K.npy", corrfull_x)
        np.save(file_path + "corr_before_mean_Y_" + str(T) + "_K.npy", corrfull_y)
        np.save(file_path + "corr_before_mean_Z_" + str(T) + "_K.npy", corrfull_z)
        
        # comb = combinations(np.arange(num_seeds), num_seeds - 1)
        
        k_all, k_all_x, k_all_y, k_all_z = [], [], [], []
        corr_all, corr_all_x, corr_all_y, corr_all_z = [], [], [], []

        for i in range(corrfull_x.shape[1]):
            print("list:", i)
            corrfull_x_test, corrfull_y_test, corrfull_z_test = corrfull_x[:, i], corrfull_y[:, i], corrfull_z[:, i]
            print("corrfull_x_test:", corrfull_x_test.shape)
            
            r2_x = corrfull_x_test #np.mean(corrfull_x_test, axis=1)
            r2_y = corrfull_y_test #np.mean(corrfull_y_test, axis=1)
            r2_z = corrfull_z_test #np.mean(corrfull_z_test, axis=1)
            r2 = (r2_x + r2_y + r2_z)/3
            
            print("r2_x:", r2_x.shape)

            corr_all.append(r2)
            corr_all_x.append(r2_x)
            corr_all_y.append(r2_y)
            corr_all_z.append(r2_z)
            
            intg_x2 = (r2_x[:-1] + r2_x[1:])*(sample_freq/2)
            intg_y2 = (r2_y[:-1] + r2_y[1:])*(sample_freq/2)
            intg_z2 = (r2_z[:-1] + r2_z[1:])*(sample_freq/2)
            intg2 = (r2[:-1] + r2[1:])*(sample_freq/2)
            
            print("intg_x2:", intg_x2.shape)
            
            CSum_x2 = np.cumsum(intg_x2, axis=0) 
            CSum_y2 = np.cumsum(intg_y2, axis=0) 
            CSum_z2 = np.cumsum(intg_z2, axis=0) 
            CSum2 = np.cumsum(intg2, axis=0) 
            
            print("CSum_x2:", CSum_x2.shape)
            
            k_new_x2 = CSum_x2 * (Unit_conv/sample_freq)
            k_new_y2 = CSum_y2 * (Unit_conv/sample_freq)
            k_new_z2 = CSum_z2 * (Unit_conv/sample_freq)
            k_new2 = CSum2 * (Unit_conv/sample_freq)
            
            print("k_new_x2:", k_new_x2.shape)

            k_all.append(k_new2)
            k_all_x.append(k_new_x2)
            k_all_y.append(k_new_y2)
            k_all_z.append(k_new_z2)
        
        k_all = np.asarray(k_all)
        k_all_x, k_all_y, k_all_z = np.asarray(k_all_x), np.asarray(k_all_y), np.asarray(k_all_z)
        
        print("Temp:", temp, "k_all shape:", k_all.shape)

        print("Temp:", temp, "k_all_directions shape:", np.asarray(k_all_x).shape, np.asarray(k_all_y).shape, np.asarray(k_all_z).shape)
        
        corr_time = np.arange(1, k_new_x2.shape[0]+1) * (sample_freq*deltaT)
        fig1 = file_path + "k_across_seeds_" + str(T) + ".png"

        start_corr_time = int(3e-12/sample_freq/deltaT)
        end_corr_time = int(4e-12/sample_freq/deltaT)

        k_temp = (k_all_x + k_all_y + k_all_z)/3
        k_mean = np.mean(k_temp[:, start_corr_time:end_corr_time], axis=1)
        kx_mean = np.mean(k_all_x[:, start_corr_time:end_corr_time], axis=1)
        ky_mean = np.mean(k_all_y[:, start_corr_time:end_corr_time], axis=1)
        kz_mean = np.mean(k_all_z[:, start_corr_time:end_corr_time], axis=1)

        kmean_val = calculate_combination_means(k_mean)
        kx_val = calculate_combination_means(kx_mean)
        ky_val = calculate_combination_means(ky_mean)
        kz_val = calculate_combination_means(kz_mean)

        print(sdir, "Temp:", temp, "mean k each seed:", np.round(k_mean, 3), np.round(np.mean(kmean_val), 3), np.round(np.std(kmean_val), 3))
        print(sdir, "Temp:", temp, "directional k values:", np.round(np.mean(kx_val), 3), np.round(np.std(kx_val), 3), np.round(np.mean(ky_val), 3), np.round(np.std(ky_val), 3), np.round(np.mean(kz_val), 3), np.round(np.std(kz_val), 3))

        plt.figure()
        for seed_ind in range(len(seed_range)):
            plt.plot(corr_time, k_all[seed_ind, :], label='seed ' + str(seed_range[seed_ind]))
        plt.plot(corr_time, np.mean(k_all, axis=0), 'k', label='avg')
        plt.legend()
        plt.ylabel("k across seeds(W/mK)")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig1, bbox_inches='tight')
        # plt.show()
        
        fig2 = file_path + "k_across_direction_X_" + str(T) + ".png"
        plt.figure()
        for seed_ind in range(len(seed_range)):
            plt.plot(corr_time, k_all_x[seed_ind, :])
        # plt.plot(corr_time, np.mean(np.asarray(k_all_y)[:, :], axis=0), label="k_y")
        # plt.plot(corr_time, np.mean(np.asarray(k_all_z)[:, :], axis=0), label="k_z")
        plt.plot(corr_time, np.mean(np.asarray(k_all_x)[:, :], axis=0), 'k', label="k_x_avg")
        plt.legend()
        # plt.ylim(-50, 50)
        plt.ylabel("k across direction (W/mK)")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig2, bbox_inches='tight')
        # plt.show()

        fig2 = file_path + "k_across_direction_Y_" + str(T) + ".png"
        plt.figure()
        for seed_ind in range(len(seed_range)):
        # plt.plot(corr_time, np.mean(np.asarray(k_all_x)[:, :], axis=0), label="k_x")
            plt.plot(corr_time, k_all_y[seed_ind, :])
        # plt.plot(corr_time, np.mean(np.asarray(k_all_z)[:, :], axis=0), label="k_z")
        plt.plot(corr_time, np.mean(np.asarray(k_all_y)[:, :], axis=0), 'k', label="k_y_avg")
        # plt.ylim(-50, 50)
        plt.legend()
        plt.ylabel("k across direction (W/mK)")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig2, bbox_inches='tight')
        # plt.show()

        fig2 = file_path + "k_across_direction_Z_" + str(T) + ".png"
        plt.figure()
        for seed_ind in range(len(seed_range)):
        # plt.plot(corr_time, np.mean(np.asarray(k_all_x)[:, :], axis=0), label="k_x")
        # plt.plot(corr_time, np.mean(np.asarray(k_all_y)[:, :], axis=0), label="k_y")
            plt.plot(corr_time, k_all_z[seed_ind, :])
        plt.plot(corr_time, np.mean(np.asarray(k_all_z)[:, :], axis=0), 'k', label="k_z_avg")
        # plt.ylim(-50, 50)
        plt.legend()
        plt.ylabel("k across direction (W/mK)")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig2, bbox_inches='tight')
        # plt.show()
    
        np.save(file_path + "k_mean_X" + str(T) + "K.npy", k_all_x)
        np.save(file_path + "k_mean_Y" + str(T) + "K.npy", k_all_y)
        np.save(file_path + "k_mean_Z" + str(T) + "K.npy", k_all_z)

        corr_time1 = np.arange(0, k_new_x2.shape[0]+1) * (sample_freq*deltaT)
        fig3 = file_path + "HCACF_across_direction_" + str(T) + ".png"
        corr_all_x, corr_all_y, corr_all_z, corr_all = np.asarray(corr_all_x), np.asarray(corr_all_y), np.asarray(corr_all_z), np.asarray(corr_all)
        corr_all_x_norm, corr_all_y_norm, corr_all_z_norm = corr_all_x/corr_all_x[:, 0].reshape(-1, 1), corr_all_y/corr_all_y[:, 0].reshape(-1, 1), corr_all_z/corr_all_z[:, 0].reshape(-1, 1)
        corr_all_norm = corr_all/corr_all[:, 0].reshape(-1, 1)

        plt.figure()
        plt.plot(corr_time1, np.mean(corr_all_x, axis=0), label="total_corr_x")
        plt.plot(corr_time1, np.mean(corr_all_y, axis=0), label="total_corr_y")
        plt.plot(corr_time1, np.mean(corr_all_z, axis=0), label="total_corr_z")
        plt.plot(corr_time1, np.mean(corr_all, axis=0), label="total_corr")
        plt.legend()
        plt.ylabel("New HCACF across direction")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig3, bbox_inches='tight')
        # plt.show()

        fig4 = file_path + "HCACF_across_seeds_" + str(T) + ".png"
        corr_all = np.asarray(corr_all)    
        # corr_save = corr_all_norm
        arr_name = "corr_mean" + str(T) + "K.npy"
        np.save(file_path + "corr_mean_X" + str(T) + "K.npy", corr_all_x)
        np.save(file_path + "corr_mean_Y" + str(T) + "K.npy", corr_all_y)
        np.save(file_path + "corr_mean_Z" + str(T) + "K.npy", corr_all_z)

        plt.figure()
        for seed_ind in range(len(seed_range)):
            plt.plot(corr_time1, corr_all[seed_ind, :]/corr_all[seed_ind, 0], label='seed ' + str(seed_range[seed_ind]))
        # plt.plot(corr_time1, np.mean(corr_all_norm, axis=0), label="avg")
        plt.legend()
        plt.ylabel("New HCACF across seeds")
        plt.xlabel("Correlation time (s)")
        plt.savefig(fig4, bbox_inches='tight')
        # plt.show()
        
        print("Temp:", temp, "end of check 3:", (time.time() - start_time))
