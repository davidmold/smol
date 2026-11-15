// Reaction probe for Stringmol 0.2.3.4 (the version used in Hickinbotham, Stepney & Hogeweg
// 2021), spatial code path: get_sw + set_exec, then spatial_exec_step until the reaction ends.
// Same output format as spreact_dev.cpp and smol-react --steptrace.
//
//   spreact CONF PAIRS [maxsteps]
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include "randutil.h"
#include "params.h"
#include "alignment.h"
#include "instructions.h"
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"
#include "setupSM.h"
// the MT state lives in mt19937-2.c's statics; compiled into this unit so it can be saved/restored
#include "mt19937-2.c"
#undef M

int spatial_exec_step(stringPM *A, smsprun *run, s_ag *act, s_ag *pass);
void place_mol(s_ag *ag, smsprun *run, int x, int y);

struct mol { std::string name, seq; };
static std::vector<mol> M;
static std::string ident(const char *s){
	for(auto &m:M) if(m.seq==s) return m.name;
	return std::string("new(") + std::to_string(strlen(s)) + ")";
}
static s_ag *mk(stringPM &A, smsprun *run, const std::string &s, int x, int y){
	s_ag *a = A.make_ag('X');
	a->S = (char*)calloc(A.maxl0,1); strncpy(a->S,s.c_str(),A.maxl0-1); a->len = strlen(a->S);
	a->pp = A.spl->make_parents(NULL,NULL);
	A.update_lineage(a,'I',1,NULL,NULL,0);
	a->spp = A.spl->getspp(a,A.extit,A.maxl0);
	place_mol(a, run, x, y);
	return a;
}
static int pos(s_ag *act, s_ag *pass, char **p, int side){
	return side ? (int)(p[1] - act->S) : (int)(p[0] - pass->S);
}

int main(int argc, char *argv[]){
	if(argc < 3){ fprintf(stderr,"usage: spreact CONF PAIRS [maxsteps]\n"); return 1; }
	int maxsteps = argc > 3 ? atoi(argv[3]) : 20000;
	SMspp SP; stringPM A(&SP);
	smsprun *run = NULL;
	smspatial_init(argv[1], &A, &run, 1);
	{	FILE *fp = fopen(argv[2],"r"); char t[8], n[256], s[4096];
		while(fscanf(fp,"%7s %255s %4095s",t,n,s)==3) M.push_back({n,s});
		fclose(fp); }
	const int cx = run->gridx/2, cy = run->gridy/2;
	static unsigned long mt_save[N]; int mti_save;
	for(auto &x:M) for(auto &y:M){
		sgenrand(4357);
		s_ag *a = mk(A,run,x.seq,cx,cy), *b = mk(A,run,y.seq,cx,cy+1);
		align sw;
		float p = A.get_sw(a,b,&sw);
		A.set_exec(a,b,&sw);
		s_ag *act = (a->status==B_ACTIVE)?a:b, *pass = (act==a)?b:a;
		int ax=act->x, ay=act->y, px=pass->x, py=pass->y;
		printf("PAIR %s %s p=%.6f act=%s\n", x.name.c_str(), y.name.c_str(), p, (act==a)?x.name.c_str():y.name.c_str());
		int steps=0; const char *end="ran out of steps";
		while(steps<maxsteps){
			char op = *(act->i[act->it]);
			printf("T %d %c %d %d %d %d %d %d %d %d %d %d\n", steps, op?op:'0',
				act->it, pos(act,pass,act->i,act->it), act->rt, pos(act,pass,act->r,act->rt),
				act->wt, pos(act,pass,act->w,act->wt), act->ft, pos(act,pass,act->f,act->ft),
				(int)strlen(act->S), (int)strlen(pass->S));
			A.nexthead = NULL;
			act->next=act->prev=pass->next=pass->prev=NULL;
			// product placement on the grid draws random numbers; undo those draws (smol-react has no grid)
			if(op=='%'){ memcpy(mt_save,mt,sizeof mt); mti_save=mti; }
			spatial_exec_step(&A,run,act,pass);
			if(op=='%'){ memcpy(mt,mt_save,sizeof mt); mti=mti_save; }
			steps++;
			bool act_alive = run->grid[ax][ay]==act, pass_alive = run->grid[px][py]==pass;
			for(s_ag *q=A.nexthead;q;q=q->next){
				if(q==act || q==pass) continue;
				printf("PROD %s %s\n", ident(q->S).c_str(), q->S);
				run->grid[q->x][q->y]=NULL; run->status[q->x][q->y]=G_EMPTY;
			}
			if(!act_alive || !pass_alive){ end = !act_alive ? "active destroyed" : "passive destroyed"; break; }
			if(act->status!=B_ACTIVE){ end="unbound normally"; break; }
		}
		bool act_alive = run->grid[ax][ay]==act, pass_alive = run->grid[px][py]==pass;
		printf("END steps=%d %s act=%s pas=%s\n", steps, end, act_alive?act->S:"-", pass_alive?pass->S:"-");
		run->grid[ax][ay]=NULL; run->status[ax][ay]=G_EMPTY;
		run->grid[px][py]=NULL; run->status[px][py]=G_EMPTY;
		A.nexthead = NULL;
	}
	return 0;
}
