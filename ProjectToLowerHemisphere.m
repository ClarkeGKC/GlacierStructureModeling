function  Pt_out = ProjectToLowerHemisphere(Pt_in)

[n, ~] = size(Pt_in);

if n~=3
  error('ProjectToLowerHemisphere(): Input points should have shape P(3,:)')
end  

L_up = Pt_in(3,:)>0;

Pt_out = Pt_in;

Pt_out(:,L_up) = -Pt_out(:,L_up);
