function plot_FK(lat_epi, lon_epi, sites, lat, lon, pstn_N, pstn_E, ref_st, time, tec, tec_shifted, FK, slowness, stacked_vector, max_stack_time, time_index, synthetic, lla_sites, corrected_time)
    
    R = 6371;
    
    % create a meshgrid for the slowness
    numslow = size(slowness,2);
    
    RR = 3;
    CC = 2;

    figure(150)
    clf
    subplot(RR,CC,2)
    surfc( slowness, slowness, abs(FK),'edgecolor','none');
    az = 0;
    el = 90;
    view(az,el);
    shading interp

    axis xy;
    axis image;
    axis square;
    box on

    colormap ( 'jet' );
    %brighten(0.9);
    colorbar;
    
    set ( gca, 'TickDir', 'out' );
    xlabel ( 'Sx (sec/km)','FontSize',12);
    ylabel ( 'Sy (sec/km)','FontSize',12);
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot axes crosshairs.
    %
    plot3 ( slowness, zeros(1,numslow), ones(1,numslow), 'w' );
    plot3 ( zeros(1,numslow), slowness, ones(1,numslow), 'w' );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot max value.
    %
    [maxFK,r] = max(FK);
    [maxFK,c] = max(maxFK);
    
    se = slowness(c);
    sn = slowness(r(c));
    
    plot3 ( slowness(c), slowness(r(c)), ones(size(c)), 'xw',...
      'MarkerSize', 12,...
      'LineWidth', 2 );
  
    [az,dist] = cart2pol(slowness(c), slowness(r(c)));
    
    az = 90-az*180/pi;
    if az<0
        az = az + 360;
    end
    
    if synthetic == 1
        title ( [ 'Slowness = ' num2str(dist) ', Azimuth = ' num2str(az) ' (Synth)'],'FontSize',10);
    else
        title ( {['Slowness = ' sprintf('%.2f',   dist) ' s/km, Azimuth = ' sprintf('%.0f', az)], ...
                 ['Velocity = ' sprintf('%.2f', 1/dist*1000) ' m/s']},'FontSize',10);
    end
    set(gca,'FontSize',10);
    box on
    hold off;
  
    % show the stacked signal
    subplot(RR,CC,1)
    stacked_signal(:,:) = stacked_vector(r(c),c,:);
    if synthetic == 1
        title(['Max Amplitude at: ' num2str(max_stack_time./3600) ' (Synth)'], 'FontSize', 12); 
    else
        title(['Max Amplitude at: ' num2str(max_stack_time./3600)], 'FontSize', 12); 
    end
    
    xlabel('Time [Hours]', 'FontSize', 10)
    ylabel('Amplitude [TECU]', 'FontSize', 10)
    grid on
    hold on
    colors=lines(size(sites,1));
    for i = 1:size(tec,2)
        tt(:,:) = tec(:,i);
        plot(time(:,ref_st)./3600,tt,'color',colors(i,:));
    end
    hold off
    box on
	legend(upper(sites), 'location', 'northeastoutside')
    set(gca,'FontSize',10);
    
    subplot(RR,CC,3) 
    xlabel('Time [Hours]', 'FontSize', 10)
    ylabel('Amplitude [TECU]', 'FontSize', 10)
    grid on
    hold on
    for i = 1:size(tec,2)
        plot(time(:,ref_st)./3600,tec_shifted(:,i),'color',colors(i,:));
    end
    
    % plot a vertical line
    plot([time(time_index,ref_st)/3600 time(time_index,ref_st)/3600], ylim, 'r');
    hold off
    box on
    set(gca,'FontSize',10);
    axis tight

    subplot(RR,CC,5)
    hold on
    plot(time(:,ref_st)./3600,stacked_signal, 'g');
    plot(corrected_time./3600,stacked_signal, 'k');
    
    % run an FFT on the signal to obtain frequency
    % make a uniformly sampled dataset
    tu = (min(corrected_time):5:max(corrected_time))';
    su = interp1(corrected_time, stacked_signal, tu);

    [freq_response,freq_index] = freqz(su,1,size(su,1),1/5);

    pM = max(abs(freq_response)); %magnitude
    freq = freq_index(abs(freq_response)==pM); %frequency

    %[~, imax] = max(stacked_signal(corrected_time < corrected_time(time_index - 50)));
    %freq = 1./(corrected_time(time_index) - corrected_time(imax));
    title(['Frequency (Doppler corrected): ' sprintf('%.5f', freq*1000) ' mHz'])
    axis tight
    distancia = distance(lat(time_index,ref_st).*pi/180,lon(time_index,ref_st).*pi/180,lat_epi.*pi/180,lon_epi.*pi/180,'radians').*R;
    
    % colocate the SAC file with the TEC data
    if synthetic == 1 && ~isempty(sac_f)
        
        distancia_SEIS = distance(seis_lla(1).*pi/180,seis_lla(2).*pi/180,lat_epi.*pi/180,lon_epi.*pi/180,'radians').*R;
    
        % the delay estimation might not be very precise
        % the delay can be within +-6 minutes
        
        % find the propagation delay from the seismic stn to the array
        delay_seis = (distancia - distancia_SEIS)./synth_params(5)+synth_params(7);
        % GPS and UTC difference
        delay_seis = delay_seis +16;
        
        sac_time = sac_t*3600+(delay(ref_st)+delay_seis);
        % resample the data
        ts_sac=timeseries(sac_f,sac_time);
        res_sac = resample(ts_sac,time(:,ref_st));
        
        plot(res_sac.time/3600,res_sac.data.*max(stacked_signal)./max(res_sac.data),'--b');
        seis_collocated = [res_sac.time/3600 res_sac.data.*max(stacked_signal)./max(res_sac.data)];
    else
        seis_collocated = [];
    end
    
    % plot a vertical line
    plot([time(time_index,ref_st)/3600 time(time_index,ref_st)/3600], ylim, 'r');
    hold off
    xlabel('Time [Hours]', 'FontSize', 10)
    ylabel('Amplitude [TECU]', 'FontSize', 10)
    grid on
    box on
    set(gca,'FontSize',10);
    
    % Plot array at time of max TEC
    subplot(RR,CC,4)
    hold on
    for i = 1:size(tec,2)
        plot(pstn_E(:,i).*1000,pstn_N(:,i).*1000, 'LineWidth', 2);
        plot(pstn_E(time_index,i).*1000,pstn_N(time_index,i).*1000,'sb','MarkerFaceColor',colors(i,:));
    end
    hold off
    box on
    set(gca,'FontSize',12);
    
    grid on;
    xlabel('m', 'FontSize', 12);
    ylabel('m', 'FontSize', 12);
    axis equal
    title(['Array Geometry: center at ' num2str(lat(time_index,ref_st)) ' ' num2str(lon(time_index,ref_st))]);
    text(pstn_E(time_index,:)*1000 + 5000,pstn_N(time_index,:)*1000,upper(sites));
    
    subplot(RR,CC,6)
    
    if lat(time_index,ref_st) < -60
        m_proj('stereographic', 'lat', -90, 'lon', 0, 'radius', 30);
    else
        m_proj('mercator', 'lat', [lat(time_index,ref_st)-15 lat(time_index,ref_st)+15], 'lon', [lon(time_index,ref_st)-15 lon(time_index,ref_st)+15]);
    end
    m_coast('patch',[204/255 255/255 51/255], 'edgecolor','k');
    m_grid('xaxislocation','top');
    hold on
    box on
%     % plot seismic station
%     [vnadx,vnady] = m_ll2xy(seis_lla(2),seis_lla(1));
%     scatter(vnadx, vnady);
%     text(vnadx+.003, vnady,'(seismic stn)')
    
    % plot array
    [x, y] = m_ll2xy(lon(time_index,:),lat(time_index,:));
    scatter(x, y, repmat(10,1,size(x,2)), colors, 'filled');
    
    [x, y] = m_ll2xy(lon(:,:),lat(:,:));
    plot(x, y, 'k');
    
    % plot epicenter
    [x, y] = m_ll2xy(lon_epi,lat_epi);
    scatter(x, y, '*r');
    
    [x, y] = m_ll2xy(lla_sites(:,2),lla_sites(:,1));
    scatter(x, y, 'xr');
    %text(x+.003, y,upper(sites))
    
    %load('/home/demian/Dropbox/Geofisica/Iono height/eclipses/2019/center_2019.mat')
    ec = readtable('/home/demian/Dropbox/OSU/Projects/Eclipses/2020/shapefiles/Demian_2020_TSE_matlab.txt');
    pp = shaperead('/home/demian/Dropbox/OSU/Projects/Eclipses/2020/shapefiles/Umbral_path.shp');
    %S = shaperead('/home/demian/Dropbox/Geofisica/Iono height/eclipses/2017/eclipse2017_shapefiles/center17.shp');
    %[ex,ey] = m_ll2xy(ec.Var7, ec.Var6);
    %plot(ex, ey, '.k')
    % plot central path
    m_plot(pp(1).X, pp(1).Y, '-r');
    %m_plot(pp(2).X, pp(2).Y, '.r');
    m_plot(pp(3).X, pp(3).Y, '-r');
    % plot the position of the eclipse at max time
    %ts = linspace(datetime('2019-07-02 18:02:30'), datetime('2019-07-02 20:43:30'), size(S.X,2)-1)';
    %ts = linspace(datetime('2017-08-21 17:00:02'), datetime('2017-08-21 19:49:59'), size(S.X,2)-1)';
    %[~, index] = min(abs((hour(ts) + minute(ts)/60 + second(ts)/3600) - time(time_index,ref_st)/3600));
    %plot(ex(index), ey(index), 'oc', 'MarkerSize', 20)
    
    % azimuth from the event
    true_az = azimuth(lat(time_index,ref_st),lon(time_index,ref_st),lat_epi,lon_epi);
    h1 = m_quiver(lon(time_index,ref_st),lat(time_index,ref_st), 10*sind(true_az), 10*cosd(true_az), 'b');
    % recovered azimuth from stacking
    h2 = m_quiver(lon(time_index,ref_st),lat(time_index,ref_st), 10*sind(az), 10*cosd(az), 'r');
    hold off
    legend([h1 h2], ['Great circle back-azimith: ' num2str(true_az) ' (' num2str(distancia) ' km)'],['F-K Filter back-azimuth: ' num2str(az)])
    box on

end