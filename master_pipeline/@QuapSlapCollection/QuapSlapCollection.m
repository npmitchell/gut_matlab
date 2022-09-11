classdef QuapSlapCollection < handle
    % Collection of instances of Quasi-Axisymmetric Pullback 
    %   for Surface Lagrangian Pullbacks
    %
    % Example usage:
    % see example_QuapSlapCollection.m

    properties
        QuapSlaps
        outputDir
        ss = [] 
    end
    
    methods
        function QSC = QuapSlapCollection(QSs, opts)
            QSC.QuapSlaps = QSs ;
            QSC.outputDir = opts.outputDir ;
        end

        function alignByFolds(QSC, options)
            
            % Get fold positions for each dataset
            QSC.getFolds() ;
            folds = QSC.folds ;
            sA = cell(Ndat, 1) ;
            sB = cell(Ndat, 1) ;
            sC = cell(Ndat, 1) ;
            ss = cell(Ndat, 1) ;

            for qq = 1:Ndat
                % get mean position of first fold
                Apos = mean(double(folds(:, 1)-1)./(nU-1)) ;
                sA{qq} = linspace(0, Apos, folds(qq, 1)) ;

                % get mean position of second fold
                Bpos = mean(double(folds(:, 2)-1)./(nU-1)) ;
                sBtmp = linspace(Apos, Bpos, folds(qq, 2)-folds(qq, 1)+1) ;
                sB{qq} = sBtmp(2:end) ;

                % get mean position of third fold
                Cpos = mean(double(folds(:, 3)-1)./(nU-1)) ;
                sCtmp = linspace(Bpos, Cpos, folds(qq, 3)-folds(qq, 2)+1) ;
                sC{qq} = sCtmp(2:end) ;

                % get mean position of third fold
                Dpos = 1 ;
                sDtmp = linspace(Cpos, Dpos, nU-folds(qq, 3)+1) ;
                sD{qq} = sDtmp(2:end) ;

                ss{qq} = [sA{qq}, sB{qq}, sC{qq}, sD{qq}] ;
            end
            QSC.ss = ss ;
        end

        function getFolds(QSC, options)
            % Load fold locations from disk
            for ii = 1:length(QSC.QuapSlaps)
                QSC.QuapSlaps{ii}.getFeatures() ;
                tmp = QSC.QuapSlaps{ii}.features ;
                ff = tmp.folds ;
                if ii == 1
                    foldU = zeros(length(QSC.QuapSlaps), size(ff, 2)) ;
                end
                for fid = 1:size(ff, 2)
                    foldU(ii, fid) = ff(tmp.fold_onset(fid), fid) ;
                end
            end
            folds = [35,56,80;
                29,57,82;
                29,47,67];
            QSC.folds = folds ;
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