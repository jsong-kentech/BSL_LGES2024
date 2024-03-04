

% Data list (sample 1)

sample =1;


if sample ==1
    folderpath = 'G:\공유 드라이브\BSL-Data\LGES\12_6cm2_soc10_EIS # Sample 1';
    filelist = {'PEIS_C09_anode_cycle_soc30.csv','PEIS_C09_anode_cycle_soc50.csv','PEIS_C09_anode_cycle_soc90.csv';...
            'PEIS_C09_cathode_cycle_soc30.csv','PEIS_C09_cathode_cycle_soc50.csv','PEIS_C09_cathode_cycle_soc90.csv';...
            'PEIS_C09_full_cycle_soc30.csv','PEIS_C09_full_cycle_soc50.csv','PEIS_C09_full_cycle_soc90.csv'};
elseif sample ==2
    folderpath = 'G:\공유 드라이브\BSL-Data\LGES\12_6cm2_soc10_EIS # Sample 2';
        filelist = {'PEIS_C11_anode_cycle_soc30.csv','PEIS_C11_anode_cycle_soc50.csv','PEIS_C11_anode_cycle_soc90.csv';...
            'PEIS_C11_cathode_cycle_soc30.csv','PEIS_C11_cathode_cycle_soc50.csv','PEIS_C11_cathode_cycle_soc90.csv';...
            'PEIS_C11_full_cycle_soc30.csv','PEIS_C11_full_cycle_soc50.csv','PEIS_C11_full_cycle_soc90.csv'};
else
    error('no matching sample')
end




% check if anode and cathode sum up to full cell

figure(1)
for i_soc = 1:size(filelist,2)

    % plot each
    for i_acf = 1:size(filelist,1)

        z_now = load([folderpath filesep filelist{i_acf,i_soc}]);

        subplot(1,3,i_soc)
        plot(z_now(:,2),-z_now(:,3)); hold on


    end
    
    z_anode_now = load([folderpath filesep filelist{1,i_soc}]);
    z_cathode_now =load([folderpath filesep filelist{2,i_soc}]);
    if length(z_anode_now(:,1))~=length(z_cathode_now(:,1))
        error('cathode and anode eis data have different sizes')
    end
    z_combined_now = z_anode_now + z_cathode_now;

    subplot(1,3,i_soc)
    plot(z_combined_now(:,2),-z_combined_now(:,3),'o')
    
end