function antColonyAlgorithm()
%   2017-02-23 by Tungpo in 834
cityNums = 20;
interation = 1000;
antNums = 12;
alpha = 1;
beta = 2;
decayRate = 0.5;
remnant = 0.5;
initiaPheromone = 100;
%initialization
cityList = rand(2,cityNums);  %city  coordinate
pheromone = zeros(cityNums,cityNums);   % store the pheromone
pheromone = pheromone + initiaPheromone/cityNums^2;
antPath = zeros(antNums, cityNums+1);% the last one stores the total distance
best = zeros(1,cityNums+1);
best(1,cityNums+1) = 10000;
for i = 1:interation
    
    %select ants' start locations randomly
    antPath = startLocation(antPath);
    
    %calculate the probability of each path
    probList = pathProb(pheromone, cityList, alpha, beta);
    
    %complete the cycle
    antPath = cycle(antPath, probList);
    
    %caculate the distace of each path
    antPath = calDis(antPath, cityList);
    
    %refresh the pheromone
    pheromone = refresh(pheromone, antPath, decayRate, remnant,initiaPheromone);
    
    %draw the shortest path
    best = drawPath(antPath, cityList,i,best);
end

    function probList = pathProb(pheromone, cityList, alpha, beta)
        cNums = size(pheromone,1);
        probList = pheromone;
        for ii = 1:cNums
            for jj = ii:cNums
                if ii == jj 
                    probList(ii,jj) = 0;
                else
                    probList(ii,jj) =( probList(ii,jj))^alpha * (1/sqrt( (cityList(1,ii)-cityList(1,jj))^2 + (cityList(2,ii)-cityList(2,jj))^2 ))^beta;
                    probList(jj,ii) = probList(ii,jj);
                end
            end
        end
        for ii = 1:cNums
            probList(ii,:) = probList(ii,:) / sum(probList(ii,:));
        end
        
    end

    function antPath = startLocation(antPath)
        [aNums,cNums] = size(antPath);
        cNums = cNums - 1;
        for ii = 1:aNums
            num = ceil(rand(1)*cNums);
            antPath(ii,1) = num;
        end
    end


    function antPath = cycle(antPath, probList)
        aNums = size(antPath,1);
        for ii = 1:aNums
            antPath(ii,:) = complete(antPath(ii,:),probList);
        end
    end

    function path = complete(path, probList)
        cNums = size(path,2);
        cNums = cNums - 1;
        for ii = 2:cNums
            probList(:,path(1,ii-1)) = 0;
            spaceList = zeros(1,cNums+1);
            for kk = 2:cNums+1
                spaceList(1,kk) = spaceList(1,kk-1) + probList(path(1, ii-1),kk-1);
            end
            %look for next city randomly
            seed = rand(1);
            for kk = 2:cNums+1
                if seed < spaceList(1,kk)
                    next = kk -1;
                    break;
                end
            end
            path(1,ii) = next;
            probList(:,next) = 0;
            % normalization 
            for kk = 1: cNums
                if sum(probList(kk,:))
                    probList(kk,:) = probList(kk,:)/sum(probList(kk,:));
                end
            end
        end
    end

    function antPath = calDis(antPath, cityList)
        [aNums,cNums] = size(antPath);
        cNums = cNums - 1;
        for ii = 1:aNums
            dis = 0;
            for jj = 2:cNums
                dis =dis + sqrt( (cityList(1,antPath(ii,jj))-cityList(1,antPath(ii,jj-1)))^2 + (cityList(2,antPath(ii,jj))-cityList(2,antPath(ii,jj-1)))^2 );
            end
            dis = dis + sqrt( (cityList(1,antPath(ii,1))-cityList(1,antPath(ii,cNums)))^2 + (cityList(2,antPath(ii,1))-cityList(2,antPath(ii,cNums)))^2 );
            antPath(ii,cNums+1) = dis;
        end
    end

    function pheromone = refresh(pheromone, antPath, decayRate, remnant,initiaPheromone)
        [aNums,cNums] = size(antPath);
        cNums = cNums - 1;
        pheromone = pheromone - decayRate * pheromone;
        for ii = 1:aNums
            for jj = 2:cNums
                pheromone(antPath(ii,jj-1),antPath(ii,jj)) = pheromone(antPath(ii,jj-1),antPath(ii,jj)) + remnant;
                pheromone(antPath(ii,jj),antPath(ii,jj-1)) = pheromone(antPath(ii,jj-1),antPath(ii,jj));
            end
            pheromone(antPath(ii,1),antPath(ii,cNums)) = pheromone(antPath(ii,1),antPath(ii,cNums))+ remnant;
            pheromone(antPath(ii,cNums),antPath(ii,1)) = pheromone(antPath(ii,1),antPath(ii,cNums));
        end
        pheromone = pheromone ./ (sum(sum(pheromone)) /100);
    end

    function best = drawPath(antPath, cityList,i,best)
        [~,cNums] = size(antPath);
        cNums = cNums - 1;
        [M,I] = min(antPath(:,cNums+1)); 
        if M < best(1,cNums+1)
            best = antPath(I,:);
        end
        for ii = 2:cNums
            plot( [ cityList(1,best(1,ii-1)) cityList(1,best(1,ii)) ],...
                  [cityList(2,best(1,ii-1)) cityList(2,best(1,ii))],...
                 'Marker','.',...
                 'color','r');
            hold on;
        end
        plot( [cityList(1,best(1,1)) cityList(1,best(1,cNums)) ],...
              [cityList(2,best(1,1)) cityList(2,best(1,cNums))],...
              'Marker','.',...
              'color','r');
        hold off;
        pause(0.05);
        disp(['当前代数：' num2str(i)  '当前距离：' num2str(best(1,cNums+1))]);
    end

end







