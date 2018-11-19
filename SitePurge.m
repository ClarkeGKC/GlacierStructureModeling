function P_out = SitePurge(P_in, PurgeList)

[LIA, LOCB] = ismember([P_in.Name_site], PurgeList);

P_out = P_in(~LIA);



  