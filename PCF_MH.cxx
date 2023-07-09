//programma per calcolare la pair correlation function di un insieme di particelle con potenziale di Lennard-Jones e statistica di Boltzmann con condizioni al contorno periodiche, utilizzando l'algoritmo di Metropolis-Hastings
#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<ctime>
#include<cstdlib>
#include<string>
#include<cstdarg>
#include<algorithm>
#include<cstdio>
using namespace std;
//costante di Boltzmann
double k_B=1.3806*pow(10,-23);
//numero di Avogadro
double avogadro=6.0221E23;
//valore vicino a infinito da assegnare alle divergenze, da scegliere molto piu alto rispetto agli ordini tipici (così da fare escludere il dato atto di moto al metropolis hastings)
double divergenza=1E15;
//massimo numero di volte in cui reiterare un loop potenzialmente infinito prima di considerarsi bloccato; tenuto basso, per uso in funzione MH();
int Nmax_inf_loop=30;
//particella
struct atomo{
	double x, y, z, vx, vy, vz;
};
struct ene{
	double Tot, Pot, Kin;
};
//minimo tra double 
double minof(int n_args, ...){ 								//cout<<"<1|"<<endl;
	va_list ap;
	va_start(ap, n_args);
	double min = va_arg(ap, double);
	for(int i = 1; i < n_args; i++) {
		double a = va_arg(ap, double);
		if(a < min) min = a;
    	}
    	va_end(ap); 									//cout<<"|1>"<<endl;
    	return (min);
}
//funzione che stampa un vettore
void vector_cout (vector<double> v){ 							//cout<<"<2|"<<endl;
	cout<<"{"<<endl;
	for (auto a : v) {
		cout<<a<<","<<endl;
	}
	cout<<"}"<<endl;								//cout<<"|2>"<<endl;
	return;
}
//funzione che stampa vettore di atomi in 3d
void gas_print(vector<atomo> v, string filename="gas.txt"){ 				//cout<<"<3|"<<endl;
	ofstream OP(filename);
	for(auto a : v){
		OP<<a.x<<" "<<a.y<<" "<<a.z<<endl;
	}										//cout<<"|3>"<<endl;
	return;
}
//distanza (con periodic boundary conditions)
double dist(atomo a1, atomo a2, double L){  						//cout<<"<4|"<<endl;
	double distanza=sqrt(pow(minof(2, abs(a1.x-a2.x),abs(L-a1.x+a2.x)),2)+pow(minof(2, abs(a1.y-a2.y), abs(L-a1.y+a2.y)),2)+pow(minof(2, abs(a1.z-a2.z), abs(L-a1.z+a2.z)),2)); 					//cout<<"|4>"<<endl;
	//cout<<"dist="<<distanza<<endl;
	return(distanza);
}
//funzione da coordinate q(x1, y1, z1, x2, y2,..., zN, vx1, vy1,...., vyN, vzN) a vettore di struct atomo in periodic boundary condition (che vengono applicate sulla base dell'L)
vector<atomo> q2atoms(vector<double> q, double L){ 					//cout<<"<5|"<<endl;
	int size=q.size();
	vector<atomo> v;
	atomo a;
	for(int i=0; i<size/2; i++){
		if(i%3==0) {
			a.x=q[i]-floor(q[i]/L)*L;
			a.vx=q[i+size/2];
		}
		else if(i%3==1) {
			a.y=q[i]-floor(q[i]/L)*L;
			a.vy=q[i+size/2];
		}
		else {
			a.z=q[i]-floor(q[i]/L)*L;
			a.vz=q[i+size/2];
			v.push_back(a);
		}
	}									//cout<<"|5>"<<endl;
	return(v);
}
//funzione che applica boundary condition su cubo di lato L a componente
void BC(double &pos, double L){						//cout<<"|5.1>"<<endl;
	pos-=floor(pos/L)*L;							//cout<<"|5.1>"<<endl;
	return;
}
//Potenziale di Lennard-Jones
double U_LJ(atomo a1, atomo a2, vector<double> para){ 			//cout<<"<6|"<<endl;
	double distanza=dist(a1,a2, para[4]);							//cout<<distanza<<endl;
	if(distanza!=0) {
		/*double my_return= (4*para[2]*(pow(para[3]/distanza, 12)-pow(para[3]/distanza, 6)));
		if (!isinf(my_return)) return (my_return); 
		else return divergenza; */
		return (4*para[2]*(pow(para[3]/distanza, 12)-pow(para[3]/distanza, 6)));
	}
	else {
		return divergenza; cout<<"div"<<endl;
	} 
}
//variabile globale vettore di energie totali, per aumentare l'efficienza del programma.
vector<double> Etot_vec;
vector<double> Ekin_vec;
double Etot;
double Ekin;
//energia ("K"=cinetica, "U"=potenziale, "T"=totale)
ene Energia(vector<atomo> gas, vector<double> para, double (*U)(atomo, atomo, vector<double>)){//opzione per capire quale energia restituire 									
												//cout<<"<6.1|"<<endl;
	ene myene;
	myene.Kin=0;
	myene.Pot=0;
	myene.Tot=0;
	for(int i=0;i<gas.size();i++){
		myene.Kin+=para[0]/2*(pow(gas[i].vx,2)+pow(gas[i].vy,2)+pow(gas[i].vz,2));
		for(int j=0;j<gas.size();j++){
			if(i!=j) myene.Pot+=U(gas[i],gas[j], para);
		}
	}
	myene.Pot/=2; 
	myene.Tot=myene.Pot+myene.Kin;	 	//cout<<Etot<<endl;		
									//cout<<"|6.1>"<<endl;
	Etot=myene.Tot;
	Ekin=myene.Kin;
	return myene;	
}
//statistica di Boltzmann 
long double Bo(vector<atomo> q, vector<double> para){ //para={m,T, eps, sigma, L}
	if(para[1]!=0) {									
		return (exp(-Energia(q, para, U_LJ).Tot/(k_B*para[1]))); //SOSTITUIRE U_LJ CON UN ALTRO POT. SE SI VUOLE
	}
	else return 0;
}
//struttura pdf con argomenti
struct statistica{
	long double (*pdf)(vector<atomo> q, vector<double>);
	vector<double> para;
	//punto di max e delta larghezza per metropolis (indicativamente opportuno dia p=1/3 di estrarre un punto)
	vector<atomo> q0; //valori di partenza delle coordinate
	double del_x; //delta delle distribuzioni uniformi di MH
	double del_v;
};

//generatore di vettore di valori casuali secondo una statistica s con Metropolis, N è il numero di vettori da generare
vector<vector<atomo>> MH(statistica s, int N, int burn_in=50){ //il numero di burn_in sono i punti iniziali che scarta poichè potrebbero avere distribuzione differente dall'attesa
	int good=0, bad=0;									//cout<<"<8|"<<endl;
	vector<vector<atomo>> success;
	success.push_back(s.q0);
	int dim=s.q0.size();			cout<<"\"prob\" q0: "<<s.pdf(s.q0, s.para)<<endl;
	int count=1;
	long double p_prec=s.pdf(s.q0, s.para);
	bool was_bad=false;//discerne se lo step precedente � stato di good++ o bad++, per aggiornare correttamente i vettori di energia
	while(count<N+burn_in){
		//riempimento vettori energia cinetica e totale per stime successive
		if(count>=burn_in){ //cout<<Ekin<<" "<<Ekin_vec.size()<<" "<<Etot<<" "<<Etot_vec.size()<<" "<<was_bad<<endl;
			if(was_bad==false || Ekin_vec.size()==0){
				Ekin_vec.push_back(Ekin);
				Etot_vec.push_back(Etot);
			}
			else {
				Ekin_vec.push_back(Ekin_vec[Ekin_vec.size()-1]);
				Etot_vec.push_back(Etot_vec[Etot_vec.size()-1]);
			}
		}
		was_bad=false;
		//generazione del q trial 
		int N_inf_loop=0;
		vector<atomo> q;
		long double p_trial;
		do{
			if (N_inf_loop==Nmax_inf_loop) {
					cout<<"bloccato in regione di inf/NaN!"<<endl; 
					goto continuazione;
			}
			q=success[count-1];
			int at_da_ca=rand()%dim;
			q[at_da_ca].x+=((double)rand()/RAND_MAX*2-1)*s.del_x;	BC(q[at_da_ca].x, s.para[4]); //applicazione delle Boundary conditions
			q[at_da_ca].y+=((double)rand()/RAND_MAX*2-1)*s.del_x;	BC(q[at_da_ca].y, s.para[4]);	
			q[at_da_ca].z+=((double)rand()/RAND_MAX*2-1)*s.del_x;	BC(q[at_da_ca].z, s.para[4]);
			q[at_da_ca].vx+=((double)rand()/RAND_MAX*2-1)*s.del_v;	
			q[at_da_ca].vy+=((double)rand()/RAND_MAX*2-1)*s.del_v;	
			q[at_da_ca].vz+=((double)rand()/RAND_MAX*2-1)*s.del_v;	
			p_trial=s.pdf(q, s.para); 
			N_inf_loop++;	
		} while (!isfinite(p_trial));
		//selezione di Metropolis
		if (p_trial>p_prec ) { 
			success.push_back(q); //gas_print(q2atoms(x, s.para[4]), "gas"+to_string(count)+".txt"); //vector_cout(x);
			p_prec=p_trial;
			if(count>=burn_in) good++;
		} 
		else {
			double r=(double)rand()/RAND_MAX;
			if(log(r)<log(p_trial)-log(p_prec)) { //uso logaritmo per gestire meglio valori al limite del range disponibile
				success.push_back(q);
				p_prec=p_trial; 
				if(count>=burn_in) good++;
			}
			else {
				success.push_back(success[count-1]);  //vector_cout(x);
				was_bad=true;
				if(count>=burn_in) bad++; //cout<<p_trial<<endl;
			}
		}
		count++; 				//Ek+=Energia(success[count-1], s.para, U_LJ, "K");
	}
	continuazione:
	success.erase(success.begin(), success.begin()+burn_in); 
	cout<<"x_trial accettati: "<<good<<endl<<"x_trial scartati: "<<bad<<endl;		//cout<<"|8>"<<endl;
	return(success);	
}
//Pair Correlation Function
vector<double> PCF(vector<vector<atomo>> ev, double L, double dr){				//cout<<"<9|"<<endl;
	vector<double>funzione(int(L/(2*dr)), 0);
	//vector<double>funzione_conta(int(minof(3, L)/(2*dr)), 0);
	int na=ev[0].size();
	for(int k=0; k<ev.size(); k++){
		for (int i=0; i<na; i++){
			for(int j=i+1; j<na; j++){
				double distanza=dist(ev[k][i], ev[k][j], L);
				if(int(distanza/dr)<funzione.size() && distanza>0) {
					funzione[int(distanza/dr)]+=2*pow(L,3)/(ev.size()*4*M_PI*pow(distanza,2)*na*na*dr);
					//funzione_conta[int(distanza/dr)]++;	
				}
			}
		}
	}												//cout<<"|9>"<<endl;
	//for(int i=0;i<funzione_conta.size();i++) cout<<i*dr+dr/2<<" "<<funzione_conta[i]<<endl;
	return(funzione);
}
//funzione che applica una generica funzione (dotata eventualmente di un parametro anche vettoriale) a tutti gli elementi di un vettore generico, restituendo i risultati in un vettore (ovvero applica una funzione element-wise a un vettore)
template<typename in_type, typename out_type, typename para_type>
vector<out_type> f_el_wise(vector<in_type> v_in, out_type (*f)(in_type, para_type), para_type para){  //cout<<"<10|"<<endl;
	vector<out_type> v_out;
	for(int i=0; i<v_in.size(); i++){
		v_out.push_back(f(v_in[i], para)); 
	} 												//cout<<"|10>"<<endl;
	return(v_out);
}
//funzione che stampa funzione 
void f_print(double (*f)(double), double a, double b, int N, string filename){         			//cout<<"<11|"<<endl;
	ofstream OP(filename);
	for (int i=0; i<N; i++){
		OP<<i*(double)(b-a)/N<<" "<<f(i*(double)(b-a)/N)<<endl;
	}												//cout<<"|11>"<<endl;
	return;
}
//funzione che istogramma vettore di valori
vector<vector<double>> isto(vector<double> v, int n_bins){						//cout<<"<11.1|"<<endl;
	double min=*min_element(v.begin(), v.end());	
	double max=*max_element(v.begin(), v.end());
	int dim=n_bins;
	double delta=(max-min)/dim;
	vector<vector<double>> my_isto(dim+1, vector<double>(2, 0));
	for (int i=0; i<=dim; i++){
		my_isto[i][0]=i*delta;
	}
	for(auto a : v){
		if(int((a-min)/delta)<=dim) my_isto[int((a-min)/delta)][1]++;
	}												//cout<<"|11.1>"<<endl;
	return(my_isto);										
}
vector<double> Smooth(vector<double> v, int k){							//cout<<"<11.2|"<<endl;
	int size=v.size();
	vector<double> smooth;
	for(int i=0;i<size-k;i++){
		double sum=0;
		for(int j=0; j<k; j++){
			sum+=v[i+j];
		}
		smooth.push_back(sum/k);
	}												//cout<<"|11.2>"<<endl;
	return(smooth);
}
//----------------------------------------------------------------------------------------------------------------
int main(){											//cout<<"<12|"<<endl;
//IMPOSTAZIONI INIZIALI__________________________________________________________________________________________
	srand(time(NULL));
	//numero atomi e numero configurazioni MH
	int na=pow(7, 3);//300; //dev'essere un cubo perfetto solo per il caso di posizionamento iniziale in forma di cristallo
	int N_MH=100000;
	//temperatura pressione e massa particelle
	double m=1.66E-27*55.85;
	double T=10;
	double P=1.01E05; //modificare in base a esigenze
	//parametri Lennard-Jones
	double eps=0.2*1.6E-19;//1E-20;//
	double sigma=2.4E-10;//1E-10;//
	//tipo di inizializzazione
	string tipo_volume= "gi"; //"st" produce reticolo con distanza media circa sigma, "gi" applica la legge dei gas ideali, "man" imposta volume manuale
	string conf_iniz= "cr"; //"rnd" produce inizializzazione casuale, "cr" reticolo cristallino cubico
	//identificatore per stampa (aggiunge la stringa ai titoli)
	string identif="_10K";
//...............................................................................................................
	//dimensioni cella spaziale (condizioni al contorno periodiche)
	double L_gi=pow(na*k_B*T/P, (double)1./3); //legge dei gas
	double L_st=pow(2, 1./6)*sigma*pow(na, 1./3)*2; //stabile (distanza tra adiacenti in solido crist. cubico è il min del pot per un fattore piccolo alla fine)
	double L_man=1E-09;
	double L; if(tipo_volume == "st") L=L_st; if(tipo_volume=="gi") L=L_gi; if(tipo_volume=="man") L=L_man;
	cout<<tipo_volume<<" "<<conf_iniz<<endl;
	cout<<"sigma/L: "<<sigma/L<<endl;						 	//cout<<"<12.0|"<<endl;
	cout<<"L: "<<L<<" m"<<endl;
	if(tipo_volume=="st") cout<<"potenziale tra due adiacenti in L_st: "<<4*eps*(pow(sigma/(L/pow(na, 1./3)), 12)-pow(sigma/(L/pow(na,1./3)), 6))<<endl<<"Pressione (l. g. i.): "<<na*k_B*T/pow(L,3)<<endl;
	//grana calcolo (e stampa) pcf
	double dr=sigma/20;
	cout<<"dr "<<dr<<endl;
	
	
													//cout<<"|12.1|"<<endl;
//INIZIALIZZAZIONE__________________________________________________________________________________________________
	//inizializzo la statistica, di boltzmann
	statistica bol;
	bol.pdf=Bo; //pdf
	bol.para={m,T, eps, sigma, L};
	//MODIFICO IL SEGUENTE A INIZIALIZZ DIRETTAMENTE A VEC<ATOMO>
	//sqrt(8*T*k_B*na/(M_PI*m)); //inizializzo il vettore x0 per metropolis, con posizioni e velocità
	vector<atomo> my_q0;
	if(conf_iniz=="rnd"){
		for(int i=0; i<na; i++) {
			atomo a;
			a.x=(double)rand()/RAND_MAX*L;
			a.y=(double)rand()/RAND_MAX*L;
			a.z=(double)rand()/RAND_MAX*L;
			a.vx=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
			a.vy=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
			a.vz=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
			my_q0.push_back(a);
		}
	}
	else if(conf_iniz=="cr"){
		if(pow(int(cbrt(na)),3)!=na) cout<<"na non è un cubo perfetto"<<endl;
		int counter=0;
		while(counter<na){
			int n_L=int(cbrt(na)); cout<<"nL^3: "<<pow(n_L, 3)<<endl;
			double space=L/n_L;
			for (int i=0; i<n_L; i++){
				for(int j=0; j<n_L; j++){
					for(int k=0; k<n_L; k++){
						atomo b;
						b.x=i*space;
						b.y=j*space;
						b.z=k*space;
						b.vx=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
						b.vy=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
						b.vz=((rand()%2)*2-1)*sqrt(2*k_B*T/m)/sqrt(3);
						my_q0.push_back(b);
						counter++;
					}
				}
			}
		}
	}
	else cout<<"scegliere configurazione iniziale"<<endl;
	//stringa di comodo
	string nomefile;
	//stampa configurazione iniziale
	nomefile= "pcf_configurazione_i"+identif+".txt";
	ofstream OP_conf(nomefile);
	for(int i=0; i<my_q0.size(); i++){
		OP_conf<<my_q0[i].x<<" "<<my_q0[i].y<<" "<<my_q0[i].z<<endl;
	}
	bol.q0=my_q0; //gas_print(q2atoms(my_x0, L));
	vector<double> my_del; //inizializzo i del per metropolis
	double my_del_x=0.01*L+0*sigma;//L/cbrt(na);//(W*sigma+(1-W)*L);//*** ATTENZIONE MODIFICARE *** e **** per avere una miglior convergenza di MH
	double my_del_v=sqrt(2*k_B*T/m)/sqrt(3)*2;//****
	bol.del_x=my_del_x;
	bol.del_v=my_del_v;
													//cout<<"|12.2|"<<endl;
//ESECUZIONE__________________________________________________________________________________________________________
	//applico metropolis hastings
	vector<vector<atomo>> MH_distro=MH(bol, N_MH);

	//calcolo la pcf
	vector<double> my_pcf=PCF(MH_distro, L, dr);
//ULTERIORI MISURE E STAMPA RISULTATI___________________________________________________________________________________

	//stampo la pcf
	nomefile="pcf"+identif+".txt";
	ofstream OP1(nomefile);									//cout<<"|12.2.1|"<<endl;
	//cout<<my_pcf.size()<<endl;
	for (int i=0; i < my_pcf.size();i++){
		OP1<<i*dr+dr/2<<" "<<my_pcf[i]<<endl;
	} 
	//stampo versione smoothed della pcf
	nomefile="pcf_smooth"+identif+".txt";
	ofstream OP1_s(nomefile);
	int smoothness=10;
	//cout<<Smooth(my_pcf).size()<<endl;
	for (int i=0; i < Smooth(my_pcf, smoothness).size();i++){				//cout<<"|12.3|"<<endl;
		OP1_s<<i*dr+dr*smoothness/2<<" "<<Smooth(my_pcf, smoothness)[i]<<endl;
	} 
	//calcolo e istogrammo energie e temperatura media					  cout<<"|12>"<<endl;
/*	vector<double> Etot_vec;
	double Ek=0;
	for(auto q : MH_distro){
		ene qene=Energia(q, bol.para, U_LJ);
		Etot_vec.push_back(qene.Tot);
		Ek+=qene.Kin;
	}
*/
	double Ek_sum=0;
	for (auto c : Ekin_vec){
		Ek_sum+=c;
	}
	cout<<"temperatura media: "<<2./((3*na-3)*(N_MH-50)*k_B)*Ek_sum<<endl; //CORREGGERE SE SI AUMENTA IL BURN_IN	
	//stampo istogramma di energie, traceplot delle energie, grafico di correlazione e funzione di autocorrelazione
	vector<vector<double>>isto_E=isto(Etot_vec, int(2*cbrt(N_MH)));
	nomefile="pcf_energie"+identif+".txt";
	ofstream OP2(nomefile);
	nomefile="pcf_energie_traceplot"+identif+".txt";
	ofstream OP2s(nomefile);
	nomefile="pcf_energie_correlazione"+identif+".txt";
	ofstream OP2c(nomefile);
	nomefile="pcf_energie_ACF"+identif+".txt";
	ofstream OP2fc(nomefile); 
	for(auto a : isto_E) {
		OP2<<a[0]<<" "<<a[1]<<endl;
	}
	double ac_0;
	double E_media=0;
	for(int i=0; i<Etot_vec.size();i++){
		E_media+=Etot_vec[i];
	}
	E_media/=Etot_vec.size();
	for(int i=0; i<Etot_vec.size();i++){
		ac_0+=pow(Etot_vec[i]-E_media,2);
	}
	ac_0/=Etot_vec.size();
	for(int l=1; l<Etot_vec.size(); l++){
		OP2s<<l<<" "<<Etot_vec[l]<<endl;
		OP2c<<Etot_vec[l-1]<<" "<<Etot_vec[l]<<endl;
		double sum=0;
		for(int i=0; i<Etot_vec.size()-l;i++){
			sum+=(Etot_vec[i]-E_media)*(Etot_vec[i+l]-E_media);
		}
		OP2fc<<l<<" "<<sum/(Etot_vec.size())/ac_0<<endl;
	}											//cout<<"|12.4|"<<endl;
	//stampa configurazione finale
	cout<<na<<" "<<MH_distro.front().size()<<" "<<MH_distro.back().size()<<endl;
	nomefile="pcf_configurazione_f"+identif+".txt";		
	ofstream OP_conf_fin(nomefile);
	for(int i=0; i<MH_distro.back().size(); i++){
		OP_conf_fin<<MH_distro.back()[i].x<<" "<<MH_distro.back()[i].y<<" "<<MH_distro.back()[i].z<<endl;
	}
	return 0;

}

