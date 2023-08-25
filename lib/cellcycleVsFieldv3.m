function [ cellCyclePhase,fieldCell,cellCyclePhaseCell ] = cellcycleVsFieldv3( dataCell,field )
%[ cellCyclePhase,fieldCell ] = cellcycleVsField( dataCell,field )

cellCyclePhase=cell(1,length(dataCell));
fieldCell=cell(1,length(dataCell));
%divTime is newborn cell time
for i=1:length(dataCell.(field))
    tempVector=zeros(length(dataCell.(field){i}),1);
    
    

     divisionTimes =  dataCell.divTime{i};

    tempVector(1:divisionTimes(1)-1)=NaN;
    for j=2:length(divisionTimes)
        indexBeforeCrop=  linspace(0,1,divisionTimes(j)-divisionTimes(j-1));
        %don't know exact cell division time
        tempVector( divisionTimes(j-1):(divisionTimes(j)-1) )=indexBeforeCrop;
    end
    tempVector(divisionTimes(end):length(tempVector))=NaN;
   
    cellCyclePhase{i}=tempVector;
       
        tempVector=dataCell.(field){i};
        fieldCell{i}=tempVector;
    
end

maxLength=0;
for i=1:length(cellCyclePhase)
    
    if length(cellCyclePhase{i})>maxLength
        maxLength=length(cellCyclePhase{i});
    end
    if length(fieldCell{i})>maxLength
        maxLength=length(fieldCell{i});
    end
    
end
dataCat1=zeros(length(cellCyclePhase),maxLength);
dataCat1(:)=NaN;
dataCat2=dataCat1;

for i=1:length(cellCyclePhase)
    dataCat1(i,1:length(cellCyclePhase{i}))=cellCyclePhase{i};
    dataCat2(i,1:length(fieldCell{i}))=fieldCell{i};
end
toRemove=isnan(dataCat1) | isnan(dataCat2);
dataCat1=dataCat1(~toRemove);
dataCat2=dataCat2(~toRemove);

cellCyclePhaseCell = cellCyclePhase;

cellCyclePhase=dataCat1;
fieldCell=dataCat2;

end

