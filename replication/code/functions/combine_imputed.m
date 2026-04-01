function combine_imputed(path,dfxin,pool,xarg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function put combines fiscal years for the binding and pooled ICC
% specifications. Specifically, it reads in the imputed costs, prices,
% shares, and supermarkups, and resaves into single file.
%
% Called by:
%   - main_v7.m
% Loads results from:
%   - main_supply_bind.m
%   - main_supply_bind_x.m 
% Saved results are loaded in:
%   - results_costregs.m
%   - results_baseanalysis.m
%   - results_cfanalysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Looping over selected dfx
for dfx = dfxin
    
    % Setting specification
    df = dfx/100;
    spec = main_spec(df);
    spec.pool = pool;
    
    % Locating results
    cd(strcat(path.data2,'/',spec.ssubfolder))
    
    
    if pool==1

        if xarg==1    
            load sres_bind_x_2006 mc pnash snash smbase smdevs
            sm_06 = smbase + smdevs;
            smbase_06 = 0;
        elseif xarg==0
            load sres_bind_2006 mc pnash snash smfinal smbase0
            temp = [smfinal(1) ; smfinal(1) + smfinal(2:end)];
            sm_06 = temp(1:length(spec.city_in));
            smbase_06 = smbase0;
        end
        mc_06 = mc;
        pn_06 = pnash;
        sn_06 = snash;
        
        if xarg==1
            load sres_bind_x_2007 mc pnash snash smbase smdevs
            sm_07 = smbase + smdevs;
            smbase_07 = 0;
        elseif xarg==0
            load sres_bind_2007 mc pnash snash smfinal smbase0
            temp = [smfinal(1) ; smfinal(1) + smfinal(2:end)];
            sm_07 = temp(1:length(spec.city_in));
            smbase_07 = smbase0;
        end
        mc_07 = mc;
        pn_07 = pnash;
        sn_07 = snash;
        
        if xarg==1
            load sres_bind_x_2010 mc pnash snash smbase smdevs
            sm_10 = smbase + smdevs;
            smbase_10 = 0;
        elseif xarg==0
            load sres_bind_2010 mc pnash snash smfinal smbase0
            temp = [smfinal(1) ; smfinal(1) + smfinal(2:end)];
            sm_10 = temp(1:length(spec.city_in));
            smbase_10 = smbase0;
        end
        mc_10 = mc;
        pn_10 = pnash;
        sn_10 = snash;
        
        if xarg==1
            load sres_bind_x_2011 mc pnash snash smbase smdevs
            sm_11 = smbase + smdevs;
            smbase_11 = 0;
        elseif xarg==0
            load sres_bind_2011 mc pnash snash smfinal smbase0
            temp = [smfinal(1) ; smfinal(1) + smfinal(2:end)];
            sm_11 = temp(1:length(spec.city_in));
            smbase_11 = smbase0;
        end
        mc_11 = mc;
        pn_11 = pnash;
        sn_11 = snash;
        
        % Combining data and resaving
        mc = [mc_06 ; mc_07 ; mc_10 ; mc_11];
        pnash = [pn_06 ; pn_07 ; pn_10 ; pn_11];
        snash = [sn_06 ; sn_07 ; sn_10 ; sn_11];
        
        smbase = [smbase_06 smbase_07 smbase_10 smbase_11];
        
        if xarg==1
            save sres_bind_x_combine mc pnash snash sm_06 sm_07 sm_10 sm_11;
        elseif xarg==0
            save sres_bind_combine mc pnash snash sm_06 sm_07 sm_10 sm_11 smbase;
        end
        
    end
    
    % Need this inside the loop
    cd(path.code1);
    
end


