function M_struct  = MergeSiteData(P_struct, SITE)

np = numel(P_struct);
ns = numel(SITE);

if np~=ns
  error('MergeSiteData(): SITE structure and P_struct have different number of elements (np~=ns)')
end

S_fields = fieldnames(SITE);
P_fields = fieldnames(P_struct);

for s=1:ns
  for q=1:numel(S_fields)
    eval(sprintf('M_struct(s).%s  = SITE(s).%s;', S_fields{q}, S_fields{q}))
  end  
  
  for p=1:numel(P_fields)
    eval(sprintf('M_struct(s).%s  = P_struct(s).%s;', P_fields{p}, P_fields{p}))
  end
end

