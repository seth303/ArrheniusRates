h2 = findall(groot,'Type','axes');

f3 = figure;
t = tiledlayout(f3,2,2);

ax1c = copyobj(h2(1),t);
ax2c = copyobj(h2(2),t);
ax3c = copyobj(h2(3),t);
ax4c = copyobj(h2(4),t);
ax1c(1).Layout.Tile = 4;
ax2c(1).Layout.Tile = 2;
ax3c(1).Layout.Tile = 3;
ax4c(1).Layout.Tile = 1;