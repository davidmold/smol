// Reaction probe for current Stringmol dev (spatial code path).
// Force-binds every listed molecule with every other (and itself), each taking a turn as
// initiator, exactly as the spatial sim does (AgentsAlign + ReactionSetupExecution), then runs
// ReactionExecuteOpcode_Spatial until the reaction ends. Mutation off. Prints a per-step trace so
// the run can be diffed against smol-react --steptrace.
//
//   spreact CONF PAIRS [maxsteps]
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <unistd.h>
#include "default_config.h"
#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"
#include "alignment.h"
#include "agent.h"
#include "SMspp.h"
#include "rules.h"
#include "opcodes.h"
#include "stringPM.h"
#include "sm_spatial.h"
#include "setupSM.h"

struct mol { std::string name, seq; };
static std::vector<mol> M;
static std::string ident(const char *s){
	for(auto &m:M) if(m.seq==s) return m.name;
	return std::string("new(") + std::to_string(strlen(s)) + ")";
}

static s_ag *mk(Stringmol_Spatial &A, const std::string &s, int x, int y){
	char buf[4096]; strncpy(buf, s.c_str(), sizeof buf - 1); buf[sizeof buf - 1] = 0;
	s_ag *a = AgentMakeWithSequence(buf, 'X', A.agct++, A.maxl0);
	a->pp = A.spl->ParentsMake(NULL,NULL);
	A.spl->SpeciesListUpdate(a,'I',1,NULL,NULL,0,A.timestep,A.maxl0);
	a->spp = A.spl->getspp(a,A.timestep,A.maxl0);
	AgentPlaceOnGrid(a, A.grid, x, y);
	return a;
}

static int pos(s_ag *act, s_ag *pass, char **p, int side){
	return side ? (int)(p[1] - act->S) : (int)(p[0] - pass->S);
}

int main(int argc, char *argv[]){
	if(argc < 3){ fprintf(stderr,"usage: spreact CONF PAIRS [maxsteps]\n"); return 1; }
	int maxsteps = argc > 3 ? atoi(argv[3]) : 20000;
	SMspp SP; Stringmol_Spatial A(&SP);
	smsprun *run = NULL;
	A.randseed = RandomSeedInitFromFile(argv[1],0);
	A.ConfigureFromFile(argv[1], &run, 1);
	{	FILE *fp = fopen(argv[2],"r"); char t[8], n[256], s[4096];
		while(fscanf(fp,"%7s %255s %4095s",t,n,s)==3) M.push_back({n,s});
		fclose(fp); }
	char statefn[64]; sprintf(statefn,"/tmp/spreact_mt_%d.txt",(int)getpid());
	const int cx = A.grid->gridx/2, cy = A.grid->gridy/2;
	for(auto &x:M) for(auto &y:M){
		SetRNGSeed(4357);   // same stream for every pair, as smol-react --steptrace does
		s_ag *a = mk(A,x.seq,cx,cy), *b = mk(A,y.seq,cx,cy+1);
		align sw;
		float p = AgentsAlign(a,b,&sw,A.blosum,A.swlist);
		A.ReactionSetupExecution(a,b,&sw);
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
			act->next=act->prev=pass->next=pass->prev=NULL;  // as if extracted from nowhead, like the sim
			// a cleave draws random numbers only to place the product on the grid; smol-react
			// has no grid, so undo those draws to keep the two streams in step
			if(op=='%'){ FILE *sf=fopen(statefn,"w"); MersenneTwisterPrintStatusToFile(sf); fclose(sf); }
			A.ReactionExecuteOpcode_Spatial(act,pass);
			if(op=='%') MersenneTwisterLoadState(statefn);
			steps++;
			bool act_alive = A.grid->grid[ax][ay]==act, pass_alive = A.grid->grid[px][py]==pass;
			for(s_ag *q=A.nexthead;q;q=q->next){
				if(q==act || q==pass) continue;
				printf("PROD %s %s\n", ident(q->S).c_str(), q->S);
				A.grid->grid[q->x][q->y]=NULL; A.grid->status[q->x][q->y]=G_EMPTY;
			}
			if(!act_alive || !pass_alive){ end = !act_alive ? "active destroyed" : "passive destroyed"; break; }
			if(act->status!=B_ACTIVE){ end="unbound normally"; break; }
		}
		bool act_alive = A.grid->grid[ax][ay]==act, pass_alive = A.grid->grid[px][py]==pass;
		printf("END steps=%d %s act=%s pas=%s\n", steps, end, act_alive?act->S:"-", pass_alive?pass->S:"-");
		A.grid->grid[ax][ay]=NULL; A.grid->status[ax][ay]=G_EMPTY;
		A.grid->grid[px][py]=NULL; A.grid->status[px][py]=G_EMPTY;
		A.nexthead = NULL;
	}
	return 0;
}
