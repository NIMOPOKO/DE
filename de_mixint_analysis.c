#include <stdio.h>
#include <stdlib.h>  // 標準ライブラリ

#define TITLE "bbob-mixint-de-sphere-lamarckian-20D"

int main(void){
    FILE *gp;

    gp = popen("gnuplot -persist","w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output './analysis/%s-divercity.png'\n",TITLE);
    // --- gnuplotにコマンドを送る --- //
	fprintf(gp, "set xrange [0.0:4.0]\n"); // 範囲の指定(省略可)
    //fprintf(gp, "set yrange [*:*] reverse\n"); // 範囲の指定(省略可)
    	
    fprintf(gp, "set xlabel 'log10(# f-evals/dimension)'\n");
    // fprintf(gp, "set ylabel 'divercity'\n");
    fprintf(gp, "set ylabel 'log10(score/max_score)'\n");
    fprintf(gp, "set grid\n");  // Enable grid lines

    //fprintf(gp, "set logscale x\n"); // x軸を対数目盛に設定 (底が10)
    fprintf(gp, "set logscale y\n"); // y2軸を対数目盛に設定 (底が10)

    // 目盛りを科学的記法で表示
    // fprintf(gp, "set format x '10^{%%L}'\n"); // y2軸の目盛りを科学的記法で表示
    fprintf(gp, "set format y '10^{%%L}'\n"); // y軸の目盛りを科学的記法で表示
    // fprintf(gp, "set format y2 '10^{%%L}'\n"); // y2軸の目盛りを科学的記法で表示


    fprintf(gp, "set xtics font ',10'\n"); // 目盛りのフォントサイズを14に設定
    fprintf(gp, "set ytics font ',10'\n"); // 目盛りのフォントサイズを14に設定
    // fprintf(gp, "set y2tics\n");
    // fprintf(gp, "set ytics nomirror\n");

    fprintf(gp,"set title '%s' font '14'\n",TITLE);
    fprintf(gp,"set key right top\n");  // Legend positioning
    fprintf(gp,"plot ");
    fprintf(gp,"'./data/%s.dat' using 1:2 w l lw 2 lc rgb 'blue' title 'sphere'",TITLE);
    fprintf(gp, "\n");
    fflush(gp); // バッファに格納されているデータを吐き出す（必須）
	pclose(gp);
    return 0;
}