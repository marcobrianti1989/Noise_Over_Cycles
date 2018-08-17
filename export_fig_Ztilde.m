function export_fig_Ztilde(export_figure1)

if export_figure1 == 1
      % Create the correct path
      base_path = pwd;
      warning off
      if exist([base_path '\Figures'], 'dir')
            cd([base_path '\Figures']) %for Microsoft
      else
            cd([base_path '/Figures']) %for Mac
      end
      if exist([base_path '\Export_Fig'], 'dir')
            addpath([base_path '\Export_Fig']) %for Microsoft
      else
            addpath([base_path '/Export_Fig']) %for Mac
      end
      warning on
      export_fig(['Noise_Shocks_SPF_GDPgrowth_revisions.pdf'])
      close 
      cd(base_path) %back to the original path
end

end