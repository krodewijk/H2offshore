props;

g.resetMask();

clust = Cluster();
clust.populate(g, 300, 200);
%clust.connect2shore(g, property.shoreConnectionPoint);
clust.plot(g);