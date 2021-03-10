props;

g.resetMask();

clust = Cluster();
clust.populate(g, 600, 300);
clust.connect2backbone(g,50);
%clust.connect2shore(g, property.shoreConnectionPoint);

clust2 = Cluster();
clust2.populate(g, 700, 300);
clust2.connect2backbone(g,50);

%%
clust.plot(g);
clust2.plot(g)
g.plot_pipes();