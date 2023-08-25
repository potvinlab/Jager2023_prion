function [dataCell,outputArg2] = prionPartErrorSignTot_keepOutliers(dataCell,lossIndex,lossTime, neverIndex,autofluo)

bleedfactor = 0;

% start with cells that never had the prion
motherYFPcat = [];
sisterYFPcat = [];
motherAreacat = [];
sisterAreacat = [];
frameCat = [];
dataCell.Q = dataCell.frameNo;


% compatibility with prog vs sister

if ~isfield(dataCell, 'prog')
   dataCell.prog.frameNo = dataCell.frameNo;
   dataCell.prog.YFPTot = dataCell.sisterYFPTot;
   dataCell.prog.Area = dataCell.sisterArea;
end


for i=1:length(neverIndex)
    cellIndex = neverIndex(i);
    
    dataCell.Q{cellIndex}(:) = NaN;
    
    
   for j=1:length(dataCell.divTime{cellIndex})
       

      divisionAtFrame = dataCell.divTime{cellIndex}(j);
      positionInMother = find(dataCell.frameNo{cellIndex}==divisionAtFrame,1,'first');
      positionInSister = find(dataCell.prog.frameNo{cellIndex}==divisionAtFrame,1,'first');
      
      
      
      if  ~isempty(positionInSister)% && divisionAtFrame <15
         
          
         frameCat = cat(1,frameCat,divisionAtFrame); 
         motherYFPcat = cat(1,motherYFPcat,dataCell.YFPTot{cellIndex}(positionInMother)-autofluo);
         motherAreacat = cat(1,motherAreacat,dataCell.Area{cellIndex}(positionInMother));
         if positionInSister > length(dataCell.prog.YFPTot{cellIndex})
            cellIndex
         end

         
         sisterYFPcat = cat(1,sisterYFPcat,dataCell.prog.YFPTot{cellIndex}(positionInSister)-autofluo - bleedfactor*motherYFPcat(end));
         sisterAreacat = cat(1,sisterAreacat,dataCell.prog.Area{cellIndex}(positionInSister));
          
         if abs(motherAreacat(end)-sisterAreacat(end))/motherAreacat(end)  < 0.25%  && motherYFPcat(end) <500*325 && sisterYFPcat(end)<500*325
         dataCell.Q{cellIndex}(positionInMother) = -((sisterYFPcat(end) - motherYFPcat(end)))/(sisterYFPcat(end) + motherYFPcat(end));
         end
      end
       
       
       
   end
    
  
    
  
end




motherYFPcat = max(0,motherYFPcat);
sisterYFPcat = max(0,sisterYFPcat);
goodIndices = abs(motherAreacat-sisterAreacat)./motherAreacat  < 0.5  & motherYFPcat <500*325 & sisterYFPcat<500*325;
avgArea = (sisterAreacat(goodIndices) + motherAreacat(goodIndices))/2;
xToPlot = (sisterYFPcat(goodIndices) + motherYFPcat(goodIndices));
yToPlot = abs((sisterYFPcat(goodIndices) - motherYFPcat(goodIndices)));



outNever.motherYFPcat = motherYFPcat;
outNever.motherAreacat = motherAreacat;
outNever.sisterYFPcat = sisterYFPcat;
outNever.sisterAreacat = sisterAreacat;

motherYFPcatBefore= [];
motherAreacatBefore = [];
sisterYFPcatBefore = [];
sisterAreacatBefore = [];
frameCatBefore = [];

motherYFPcatAfter= [];
motherAreacatAfter = [];
sisterYFPcatAfter = [];
sisterAreacatAfter = [];
frameCatAfter=[];
for i=1:length(lossIndex)
    cellIndex = lossIndex(i);
        dataCell.Q{cellIndex}(:) = NaN;

   for j=1:length(dataCell.divTime{cellIndex})

      divisionAtFrame = dataCell.divTime{cellIndex}(j);
      positionInMother = find(dataCell.frameNo{cellIndex}==divisionAtFrame,1,'first');
      positionInSister = find(dataCell.prog.frameNo{cellIndex}==divisionAtFrame,1,'first');
      
      if  ~isempty(positionInSister)
          
          if divisionAtFrame <=lossTime(i)
                       frameCatBefore = cat(1,frameCatBefore,divisionAtFrame); 

              motherYFPcatBefore = cat(1,motherYFPcatBefore,dataCell.YFPTot{cellIndex}(positionInMother)-autofluo);
              motherAreacatBefore = cat(1,motherAreacatBefore,dataCell.Area{cellIndex}(positionInMother));
              
              sisterYFPcatBefore = cat(1,sisterYFPcatBefore,dataCell.prog.YFPTot{cellIndex}(positionInSister)-autofluo- bleedfactor*motherYFPcatBefore(end));
              sisterAreacatBefore = cat(1,sisterAreacatBefore,dataCell.prog.Area{cellIndex}(positionInSister));
              
              if abs(motherAreacatBefore(end)-sisterAreacatBefore(end))/motherAreacatBefore(end)  < 0.25  && motherYFPcatBefore(end) <500*325 && sisterYFPcatBefore(end)<500*325
                  dataCell.Q{cellIndex}(positionInMother) = -((sisterYFPcatBefore(end) - motherYFPcatBefore(end)))/(sisterYFPcatBefore(end) + motherYFPcatBefore(end));
              end
              
          else
                       frameCatAfter = cat(1,frameCatAfter,divisionAtFrame); 

              motherYFPcatAfter = cat(1,motherYFPcatAfter,dataCell.YFPTot{cellIndex}(positionInMother)-autofluo);
              motherAreacatAfter = cat(1,motherAreacatAfter,dataCell.Area{cellIndex}(positionInMother));
              
              sisterYFPcatAfter = cat(1,sisterYFPcatAfter,dataCell.prog.YFPTot{cellIndex}(positionInSister)-autofluo- bleedfactor*motherYFPcatAfter(end));
              sisterAreacatAfter = cat(1,sisterAreacatAfter,dataCell.prog.Area{cellIndex}(positionInSister));
              if abs(motherAreacatAfter(end)-sisterAreacatAfter(end))/motherAreacatAfter(end)  < 0.25  && motherYFPcatAfter(end) <500*325 && sisterYFPcatAfter(end)<500*325
                  dataCell.Q{cellIndex}(positionInMother) = -((sisterYFPcatAfter(end) - motherYFPcatAfter(end)))/(sisterYFPcatAfter(end) + motherYFPcatAfter(end));
              end
          end
      end
      
      
      
   end
   
  
    
  
end
conversionFactor = 1;
motherYFPcatAfter = max(0,motherYFPcatAfter);
sisterYFPcatAfter = max(0,sisterYFPcatAfter);

motherYFPcatBefore = max(0,motherYFPcatBefore);
sisterYFPcatBefore = max(0,sisterYFPcatBefore);


goodIndices = abs(motherAreacatBefore-sisterAreacatBefore)./motherAreacatBefore  < 0.25  & motherYFPcatBefore <500*326 & sisterYFPcatBefore<500*325;
avgArea = (sisterAreacatBefore(goodIndices) + motherAreacatBefore(goodIndices))/2;
xToPlot = (sisterYFPcatBefore(goodIndices) + motherYFPcatBefore(goodIndices))/conversionFactor;
yToPlot = abs((sisterYFPcatBefore(goodIndices) - motherYFPcatBefore(goodIndices)))/conversionFactor;



goodIndices = abs(motherAreacatAfter-sisterAreacatAfter)./motherAreacatAfter  < 0.25  ;%& motherYFPcatAfter <500*325 & sisterYFPcatAfter<500*325;
avgArea = (sisterAreacatAfter(goodIndices) + motherAreacatAfter(goodIndices))/2;
xToPlot = (sisterYFPcatAfter(goodIndices) + motherYFPcatAfter(goodIndices))/conversionFactor;
yToPlot = ((sisterYFPcatAfter(goodIndices) - motherYFPcatAfter(goodIndices))).^2/conversionFactor;

end

