function IC = toeoff(KIN, side)

axis=3;

%FIND POSITIVE PEAKS
[pp_to_l,indp_to_l]=findpeaks(KIN.Pos.(side).TOE(axis,50:end),'MinPeakDistance',50);
 
 
 
%SET START AND SHIFT ACCORDING TO IT
start=indp_to_l(2); %+round((indp_to_l(2)-indp_to_l(1)));
shift=round(4*(indp_to_l(2)-indp_to_l(1))/5);
 
%FIND TO POINT USING INVERSED CURVE
[p_to_l,ind_to_l]=findpeaks(-KIN.Pos.(side).TOE(axis,start:end),'MinPeakDistance',shift);
ind_to_l=ind_to_l+start;
 
 
[pp_to_l,indp_to_l]=findpeaks(KIN.Pos.(side).TOE(axis,:),'MinPeakDistance',shift);

IC = ind_to_l;
end