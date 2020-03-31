% plot_all_cluster_size_distributions_subplots_from_xlsx.m

% STILL EDITING

% logicals
l_assemble = 1;
l_ceil = 1;
l_scale = 0; % currently scale by tot pop, also can scale by mean pop below
l_logmean = 0;  % newest version logmean not supported! only meanlog.

strains = {'a01','a02','ent','pls','psd','z20','z36','z20-che','z20-mot'};

strain_ids = [1,2,3,4,5,6,7,8,9];

tot_pop_thresh = 0*100;

% plot params
marker_cell = {'o','o','o','o','o','o','o'};
color_cell = {[0.8392 0.3373 0.2353],[0.8000 0.3137 0.4078],[0.9294 0.5176 0.2510], ...
    [0.4157 0.3608 0.6196], [0.6353 0.8000 0.2431],  [0.4 0.4 0.4], [0.1176 0.5765 0.7765], [0 1 1], [1 0 1]};
marker_size = 14;
line_width = 2;

% load data
if ~exist('T','var')
    %T = readtable('/c/Users/Brandon/Documents/Gutz/clusters/cluster_paper/data/combined_cluster_sizes.xlsx');
    T = readtable('C:\Users\Brandon\Documents\Gutz\clusters\cluster_paper\data\combined_cluster_sizes.xlsx');
end

cluster_sizes = table2array(T(:,4:end));

% fixed bins
num_bins = 10+1;
if l_scale
    bins = logspace(-2,2,num_bins);
    %bins = logspace(-4,1,num_bins);

else
    bins = logspace(0,4.2,num_bins);
end
binwidths = diff(bins);
bins = bins(1:end-1);
num_bins = num_bins - 1;

% assemble prob dens
if l_assemble
    prob_dens_cell = cell(numel(strain_ids),1);
    tot_pop_cell = cell(numel(strain_ids),1);
    num_clusters_cell = cell(numel(strain_ids,1));
    count_cell = cell(numel(strain_ids,1));
    
    for s = 1:numel(strain_ids)
          
        this_prob_dens_mat = [];
        this_tot_pop_mat = [];
        this_num_clusters_mat = [];
        this_count_mat = [];
        
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
            
            if l_scale
                these_cluster_sizes = these_cluster_sizes./mean(these_cluster_sizes);
                %these_cluster_sizes = these_cluster_sizes./sum(these_cluster_sizes);
                
            end
            
            [counts,~] = hist(these_cluster_sizes,bins);
            %counts(counts==0) = 1;
            
            this_prob_dens_mat = [this_prob_dens_mat; counts./sum(counts)./binwidths];
            this_tot_pop_mat = [this_tot_pop_mat; sum(these_cluster_sizes)];
            this_num_clusters_mat = [this_num_clusters_mat; numel(these_cluster_sizes)];
            this_count_mat = [this_count_mat; counts];
            
        end
        
    
        prob_dens_cell{s} = this_prob_dens_mat;
        tot_pop_cell{s} = this_tot_pop_mat;
        num_clusters_cell{s} = this_num_clusters_mat;
        count_cell{s} = this_count_mat;
    
    end
    
end



%% plot

figure('position',  [488.0000   41.8000  915.4000  740.8000]); hold on;
legend_cell = cell(1,numel(strain_ids));
subplot_inds = [1,2,3,4,5,7,6,8,9];
title_cell = {'\it{A. veronii}', '\it{A. caviae}', '\it{E. cloacae}', '\it{Plesiomonas}',...
    '\it{P. mendocina}', '\it{V. cholerae} (20)','\it{V. cholerae} (36)','\it{V. cholerae} (20) {\Delta}che', '\it{V. cholerae} (20) {\Delta}mot'};
integrated_diffs = [];
for s = 1:numel(strain_ids)
    
    subplot(3,3,subplot_inds(s)); hold on;
    
    this_prob_dens_mat = prob_dens_cell{s};
    this_tot_pop = tot_pop_cell{s};
    this_num_clusters = num_clusters_cell{s};
    
    % tot pop threshhold
    this_prob_dens_mat(this_tot_pop < tot_pop_thresh,:) = [];
    

if l_logmean
 %%%%%%%%%%%%%%%%%%% log(mean(x)) %%%%%%%%%%%%%%%%  
    this_std_prob_dens = zeros(1,size(this_prob_dens_mat,2));
    this_mean_prob_dens = zeros(1,size(this_prob_dens_mat,2));
    for n = 1:size(this_prob_dens_mat,2)
        not_inf_ids = find(~isinf(this_prob_dens_mat(:,n)));
        this_mean_prob_dens(n) = log10(sum(this_prob_dens_mat(not_inf_ids,n))./numel(not_inf_ids));
        this_std_prob_dens(n) = std(this_prob_dens_mat(not_inf_ids,n))./this_mean_prob_dens(n);%./sqrt(numel(not_inf_ids));
    end
    
else
 %%%%%%%%%%%%%%%%%%% mean(log(x)) %%%%%%%%%%%
    this_log_prob_dens_mat = log10(this_prob_dens_mat);
    this_std_prob_dens_bio = zeros(1,size(this_log_prob_dens_mat,2));
    %this_mean_prob_dens = zeros(1,size(this_log_prob_dens_mat,2));
    for n = 1:size(this_log_prob_dens_mat,2)
        not_inf_ids = find(~isinf(this_log_prob_dens_mat(:,n)) & ~isnan(this_log_prob_dens_mat(:,n)));
        %this_mean_prob_dens(n) = sum(this_log_prob_dens_mat(not_inf_ids,n))./numel(not_inf_ids);
        this_std_prob_dens_bio(n) = std(this_log_prob_dens_mat(not_inf_ids,n))./sqrt(numel(not_inf_ids));
    end
    
end
   
%%%%%%%%%%%%%%%%% sampling uncertainty %%%%%%%%%%%%%%%
    this_count_mat = count_cell{s};
    this_uncertainty_count_mat = sqrt(this_count_mat);
    total_counts_mat = repmat(sum(this_count_mat,2),1,size(this_count_mat,2));
    this_log_prob_dens_mat = log10(this_count_mat./total_counts_mat./repmat(binwidths,size(this_count_mat,1),1));
    % add contribution from total counts?
    this_uncertainty_log_prob_dens_mat = sqrt( (this_uncertainty_count_mat./this_count_mat).^2 + ...
        (repmat(binwidths,size(this_count_mat,1),1)./size(this_count_mat,2)./total_counts_mat./repmat(binwidths,size(this_count_mat,1),1)).^2.*repmat(sum(this_uncertainty_count_mat.^2,2),1,size(this_count_mat,2)));
    
    
    this_mean_prob_dens = zeros(1,size(this_log_prob_dens_mat,2));
    this_std_prob_dens_sampling = zeros(1,size(this_log_prob_dens_mat,2));
    for n = 1:size(this_log_prob_dens_mat,2)
        not_inf_ids = find(~isinf(this_log_prob_dens_mat(:,n)) & ~isnan(this_log_prob_dens_mat(:,n)));
        this_mean_prob_dens(n) = sum(this_log_prob_dens_mat(not_inf_ids,n)./this_uncertainty_log_prob_dens_mat(not_inf_ids,n))./sum(1./this_uncertainty_log_prob_dens_mat(not_inf_ids,n));
        this_std_prob_dens_sampling(n) = sqrt((1./numel(this_uncertainty_log_prob_dens_mat(not_inf_ids,n)).^2).*sum(this_uncertainty_log_prob_dens_mat(not_inf_ids,n).^2));
    end
    
    this_std_prob_dens = sqrt(this_std_prob_dens_bio.^2 + this_std_prob_dens_sampling.^2);
    
    this_x = bins(~isinf(this_mean_prob_dens));
    this_sigy = this_std_prob_dens(~isinf(this_mean_prob_dens));
    this_y = this_mean_prob_dens(~isinf(this_mean_prob_dens));
    %this_std_prob_dens(this_mean_prob_dens - this_std_prob_dens <= 0 ) = this_mean_prob_dens(this_mean_prob_dens - this_std_prob_dens <= 0 ).*.9;
    %errorbar(bins,this_mean_prob_dens,this_std_prob_dens./sqrt(size(this_prob_dens_mat,1)),['k' marker_cell{strain_ids(s)}],'markersize',marker_size,'markerfacecolor',color_cell{strain_ids(s)},'linewidth',line_width);
    %plot(bins,this_mean_prob_dens,'color',color_cell{strain_ids(s)},'linewidth',4)
    %shadedErrorBar(bins,this_mean_prob_dens,this_std_prob_dens./sqrt(size(this_prob_dens_mat,1)),'lineProps',{'color',color_cell{strain_ids(s)},'linewidth',2},'transparent',true);
    %shadedErrorBar(this_x,this_y,this_sigy./sqrt(size(this_prob_dens_mat,1)),'lineProps',{'color',color_cell{strain_ids(s)},'linewidth',2},'transparent',true);
    % -2 line
    %A = .001;
    A = 1;
    %xline = bins(ceil(numel(bins)/4):ceil(2*numel(bins)/4));
    xline = bins;
    yline = A.*xline.^(-2);
    yline = log10(yline) - log10(yline(1)) + this_y(1);
    %if s==1
        h = plot(xline,yline,'k--','linewidth',2);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    %end
    
    errorbar(this_x,this_y,this_sigy./sqrt(size(this_prob_dens_mat,1)),'ko','markersize',marker_size,'markerfacecolor',color_cell{strain_ids(s)},'linewidth',2);
    
    legendcell{s} = strains{strain_ids(s)};


    set(gca,'fontsize',16,'linewidth',4,'xscale','log','yscale','linear','xtick',[1e0 1e2 1e4],'xminortick','off')
    if l_scale
        %  xlabel('scaled cluster size','fontsize',24)
    else
        %   xlabel('n (number of cells)','fontsize',24)
    end
    %ylabel('log_{10}p(n)','fontsize',24)
    %legend(legendcell,'location','ne','fontsize',16);
    axis([.1 1e4 -10 0])
    %title(strains{strain_ids(s)},'fontsize',16);
    title(title_cell{strain_ids(s)},'fontsize',16);

    if s==8
        xlabel('n (number of cells)','fontsize',16)
    end
    
    if s==4
        ylabel('log_{10}p(n)','fontsize',16)
    end
    
    
    this_diff = this_mean_prob_dens - yline;
    this_diff(isinf(this_mean_prob_dens)) = [];
    this_dx = binwidths;
    this_dx(isinf(this_mean_prob_dens)) = [];
    
    this_integrated_diff = sum(this_diff.*this_dx)./max(bins);
    integrated_diffs = [integrated_diffs, this_integrated_diff];
end
