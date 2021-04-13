classdef QuapSlapCollection < handle
    % Collection of instances of Quasi-Axisymmetric Pullback 
    %   for Surface Lagrangian Pullbacks
    %

    properties
        QuapSlaps
        outputDir
    end
    
    methods
        function QSC = QuapSlapCollection(QSs, opts)
            QSC.QuapSlaps = QSs ;
            QSC.outputDir = opts.outputDir ;
        end
        
        function collectStrainRates(QSC, options)
            % Collect strain rate kymographs for all QS objects in
            % collection.
            % 
            % Parameters
            % ----------
            % QSC : class instance
            % options : struct
            
            trsK = cell(length(QSC.QuapSlaps), 1) ;
            dvsK = cell(length(QSC.QuapSlaps), 1) ;
            thsK = cell(length(QSC.QuapSlaps), 1) ;
            for qq = 1:length(QSC.QuapSlaps)
                QS = QSC.QuapSlaps{qq} ;
                t0Pathline = QS.t0set() ;
                mKPDir = sprintf(QS.dir.strainRate.pathline.root, t0Pathline) ;
                datdir = fullfile(mKPDir, 'measurements') ;
                apKymoFn = fullfile(datdir, 'apKymographPathlineStrainRate.mat') ;
                load(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
                featureIDs = QS.getPathlineFeatureIDs('vertices', struct()) ;
                folds = load(QS.fileName.fold) ;
                fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;
                trsK{qq} = 0.5*tr_apM ;
                dvsK{qq} = dv_apM ;
                thsK{qq} = th_apM ;
            end

            close all
            set(gcf, 'visible', 'off')
            imagesc((1:nU)/nU, tps, trK)
            caxis([-climits{pp}, climits{pp}])
            colormap(bbr256)

            % Plot fold identifications
            hold on;
            fons1 = max(1, fons(1)) ;
            fons2 = max(1, fons(2)) ;
            fons3 = max(1, fons(3)) ;
            t1ones = ones(size(tps(fons1:end))) ;
            t2ones = ones(size(tps(fons2:end))) ;
            t3ones = ones(size(tps(fons3:end))) ;
            tidx0 = QS.xp.tIdx(t0) ;

            % OPTION 1: use identified div(v) < 0
            plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
            plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
            plot(featureIDs(3) * t3ones / nU, tps(fons3:end))

            % Titles 
            title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            cb = colorbar() ;

            % title and save
            ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
            disp(['saving ', fn])
            export_fig(fnout, '-png', '-nocrop', '-r200')   
        end
    end
end