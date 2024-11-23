clear all
clc

tec_project = project('/home/demian/Dropbox/OSU/Projects/Eclipses/2020/Campaña/data/config.ini');

tec_project.load_files();

for i = 1:size(tec_project.rinex_files, 1)
    files = tec_project.rinex_files{i, 1};
    orbit = tec_project.rinex_files{i, 3};
    ionex = tec_project.rinex_files{i, 4};
    
    for j = 1:size(files, 1)
        
        disp([' >> Loading rinex ' files{j,1}])
        if exist(files{j, 1}, 'file')
            stn = station(files{j, 1}, tec_project.options);
            disp(' -- Estimating TEC...')
            tec_obj = tec(stn.name, stn.rinex.observables, stn.rinex.info, orbit, ionex, tec_project.options);
            tec_obj.save(files{j, 2}, files{j, 3})
        else
            disp([' -- Could not find ' files{j, 1} ' -> skipping'])
        end
    end
end