% plot_all_cluster_size_cumdist_oneaxis_from_xlsx.m

% logicals
l_assemble = 1;
l_ceil = 1;

strains = {'a01','a02','ent','pls','psd','z20','z36','z20-che','z20-mot'};

strain_ids = [1,2,3,4,5,6,7,8,9];
color_ids = [5,4,3,1,9,2,8,7];

tot_pop_thresh = 0*100;

% plot params
marker_cell = {'o','o','o','o','o','o','o'};
color_cell = {[0.8392 0.3373 0.2353],[0.8000 0.3137 0.4078],[0.9294 0.5176 0.2510], ...
    [.72 0 .72], [0.6353 0.8000 0.2431],[0.1176 0.5765 0.7765]  , [0.4 0.4 0.4], [0 .82 .82], [0.4157 0.3608 0.6196]};
marker_size = 14;
line_width = 4;
xbar_width = .1;
reds = linspace(1,0,numel(strain_ids));
blues = ones(1,numel(strain_ids));
greens = linspace(0,1,numel(strain_ids));
% load data
if ~exist('T','var')
    %T = readtable('/c/Users/Brandon/Documents/Gutz/clusters/cluster_paper/data/combined_cluster_sizes.xlsx');
    T = readtable('C:\Users\Brandon\Documents\Gutz\clusters\cluster_paper\data\combined_cluster_sizes.xlsx');
end

cluster_sizes = table2array(T(:,4:end));



% assemble prob dens
if l_assemble
    cum_dist_cell = cell(numel(strain_ids),1);
    
    for s = 1:numel(strain_ids)
          
        if strain_ids(s) == 6
            continue
        end
        
        tmp_sizes = [];
        
        these_rows = find(strcmp(T.strain,strains{strain_ids(s)})) ;
        
            
        for r = 1:numel(these_rows)
            these_cluster_sizes = cluster_sizes(these_rows(r),:);
            these_cluster_sizes(isnan(these_cluster_sizes)) = [];
            
            % z20
            if strain_ids(s)==6
                
                these_cluster_sizes = [these_cluster_sizes, ones(1,ceil(max(these_cluster_sizes)))];
                these_cluster_sizes(these_cluster_sizes==max(these_cluster_sizes)) = [];
            end
            
            if l_ceil
                these_cluster_sizes = ceil(these_cluster_sizes);
            end
            
           
            
            tmp_sizes = [tmp_sizes, these_cluster_sizes];
            
        end
        
        tmp_sizes = sort(tmp_sizes);
        this_cum_dist = zeros(1,numel(tmp_sizes));
        for n = 1:numel(this_cum_dist)
            this_cum_dist(n) = sum(tmp_sizes > tmp_sizes(n))./numel(tmp_sizes);
        end
        
        cum_dist_cell{s} = [tmp_sizes', this_cum_dist'];
        
    end
    
end



%% plot

figure; hold on;
legend_cell = cell(1,numel(strain_ids));
subplot_inds = [1,2,3,4,5,7,6,8,9];
title_cell = {'{\it{Aeromonas}} ZOR0001', '{\it{Aeromonas}} ZOR0002', '{\it{Enterobacter}} ZOR00014', '{\it{Plesiomonas}} ZOR0011',...
    '{\it{Pseudomonas}} ZWU0006', '{\it{Vibrio}} ZWU0020','{\it{Vibrio}} ZOR0036','{\it{Vibrio}} ZWU0020 {\Delta}che', '{\it{Vibrio}} ZWU0020 {\Delta}mot'};
integrated_diffs = [];
for s = 1:numel(strain_ids)
    
    if strain_ids(s) == 6
            continue
    end
        
        
    %subplot(3,3,subplot_inds(s)); hold on;
    
    these_sizes = cum_dist_cell{s}(:,1);
    these_cum_dists = cum_dist_cell{s}(:,2);
    
    
%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%% 
    thiscolor_id = find(strain_ids(s)==color_ids);
    thiscolor = [reds(thiscolor_id),greens(thiscolor_id),blues(thiscolor_id)];
    
    
    % -1 line
    xline = these_sizes;
    yline = .05*these_sizes.^(-1);
    plot(xline,yline,'k--','linewidth',3);

    % cum dist
    %plot(these_sizes,these_cum_dists,'linewidth',line_width,'color',color_cell{strain_ids(s)});
    %p = plot(these_sizes,these_cum_dists,'o','markersize',12,'markerfacecolor',color_cell{strain_ids(s)},'color',color_cell{strain_ids(s)});
    p = plot(these_sizes,these_cum_dists,'-','linewidth',4,'color',thiscolor);

    p.Color(4) = 0.5;
    %%%%%%%%%%%%%% style %%%%%%%%%%%%%
    
    set(gca,'fontsize',16,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0 1e2 1e4],'xminortick','off','yminortick','off')
    axis([1 1e4 1e-4 1e0])
    %title(title_cell{strain_ids(s)},'fontsize',16);
    
    % optional legend: would be nice to have, but difficult to get the size
    % right in matlab.
    %if s==1
    %    legendcell = {'~n^{-2}','~n^{-1}','individual fish','average'};
     %   legend(legendcell,'location','sw','fontsize',6);
    %end

    % axis labels only once
    if s==8
        xlabel('{\it{n }}(number of cells)','fontsize',16)
    end
    
    if s==4
        ylabel('{\it{P}}(size > {\it{n}})','fontsize',16)
    end
    
    %%%%%%%%%%%%%%% misc %%%%%%%%%%%%%%%
    
    
end
