#include"lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP{
namespace lifted_disjoint_paths {

LdpInstance::LdpInstance(ConfigDisjoint<>& configParameters,char delim,CompleteStructure<>* cs,size_t minTime,size_t maxTime):
		parameters(configParameters)
{

	std::cout<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile()<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile().flush();

	useTimeFrames=parameters.isRestrictFrames()||parameters.isSparsify();
	size_t maxVertex;
	if(cs==0){
		if(useTimeFrames){
			vertexGroups=VertexGroups<size_t>(configParameters,delim);
			maxVertex=vertexGroups.getMaxVertex();
		}
		else{
			maxVertex=std::numeric_limits<size_t>::max();
		}
		std::ifstream graphFile(parameters.getGraphFileName());
		readGraph(graphFile,maxVertex,delim);
		if(!parameters.isAutomaticLifted()){
			//std::vector<std::vector<bool>> secOrderDesc=automaticLifted(graph_);
			std::cout<<"Reading lifted edges from file."<<std::endl;
			std::string line;
			std::vector<std::string> strings;
			while (std::getline(graphFile, line) && !line.empty()) {
				strings = split(line, delim);
				if (strings.size() < 3) {
					throw std::runtime_error(
							std::string("Edge vertices and score expected on every edge line "));
				}

				unsigned int v = std::stoul(strings[0]);
				unsigned int w = std::stoul(strings[1]);
				if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
				//if(isReachable(v,w)&&v!=s_&&w!=t_&&v!=t_&&w!=s_){
				if(isReachable(v,w)&&v!=s_&&w!=t_&&v!=t_&&w!=s_){
					double score = std::stod(strings[2]);
					//if(secOrderDesc[v][w]){
					auto edgeTest=graphLifted_.findEdge(v,w);
					if(!edgeTest.first){
						graphLifted_.insertEdge(v, w);
						liftedEdgeScore.push_back(score);
					}
					else{
						liftedEdgeScore[edgeTest.second]=score;
					}

				}
			}

		}
		graphFile.close();
	}
	else{
		readGraphWithTime(minTime,maxTime,cs);

	}


	//desc=disjointPaths::initReachable(graph_);

	//std::vector<std::vector<bool>> secOrderDesc=automaticLifted(graph_);

	//TODO ensure that the lifted edges are added only if the reachability and non uniqueness is satisfied!

	if(parameters.isAutomaticLifted()){
		std::cout<<"Adding automatic lifted edges"<<std::endl;
		parameters.infoFile()<<"Adding automatic lifted edges"<<std::endl;
		parameters.infoFile().flush();
		for (int i = 0; i < graph_.numberOfEdges(); ++i) {
			size_t v0=graph_.vertexOfEdge(i,0);
			size_t v1=graph_.vertexOfEdge(i,1);
			if(v0!=s_&&v1!=t_){
				//	if(secOrderDesc[v0][v1]){
				graphLifted_.insertEdge(v0,v1);
				liftedEdgeScore.push_back(edgeScore[i]);

			}
		}
		std::cout<<"done"<<std::endl;
		parameters.infoFile()<<"done"<<std::endl;
		parameters.infoFile().flush();

	}


	std::cout<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
	parameters.infoFile()<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
	parameters.infoFile().flush();
	if(parameters.isSparsify()){
		sparsifyBaseGraph();
		sparsifyLiftedGraph();

	}
	else{
		//desc=initReachable(graph_,parameters);
		reachable=initReachableSet(graph_);
	}
}

void LdpInstance::readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs){

	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
	std::vector<double>& completeScore=cs->completeScore;
	VertexGroups<> vg=cs->vg;

	std::unordered_map<size_t,std::vector<size_t>> groups;

	size_t mt=minTime;
	while(vg.getGroupVertices(mt).size()==0){
		mt++;
	}
	size_t minVertex=vg.getGroupVertices(mt)[0];

	mt=maxTime-1;
	while(vg.getGroupVertices(mt).size()==0){
		mt--;
	}
	size_t maxVertex=*(vg.getGroupVertices(mt).rbegin());

	size_t numberOfVertices=maxVertex-minVertex+3;
	s_ = numberOfVertices - 2;
	t_ = numberOfVertices - 1;

	std::vector<size_t> vToGroup(numberOfVertices);
	vToGroup[s_]=0;
	vToGroup[t_]=maxTime-minTime+1;

	for (int gi = minTime; gi < maxTime; ++gi) {
		//groups[gi-minTime+1]=std::vector<size_t>();
		for(size_t v:vg.getGroupVertices(gi)){
			size_t vertex=v-minVertex;
			groups[gi-minTime+1].push_back(vertex);
			vToGroup[vertex]=gi-minTime+1;
		}
	}

	vertexGroups=VertexGroups<>(groups,vToGroup);
	vertexScore = std::vector<double>(numberOfVertices, 0);
	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
	graph_ = andres::graph::Digraph<>(numberOfVertices);

	bool useZeroInOut=false;
	for (int v = 0; v < numberOfVertices-2; ++v) {
		graph_.insertEdge(s_,v);
		if(useZeroInOut) edgeScore.push_back(0);
		else edgeScore.push_back(parameters.getInputCost());
		graph_.insertEdge(v,t_);
		if(useZeroInOut) edgeScore.push_back(0);
		else edgeScore.push_back(parameters.getOutputCost());
	}


	for (int v = minVertex; v < maxVertex; ++v) {
		for (int i = 0; i < completeGraph.numberOfEdgesFromVertex(v); ++i) {
			size_t w=completeGraph.vertexFromVertex(v,i);
			if(w>maxVertex) continue;
			size_t e=completeGraph.edgeFromVertex(v,i);
			graph_.insertEdge(v-minVertex, w-minVertex);
			edgeScore.push_back(completeScore[e]);
		}
	}
	minV=minVertex;
	maxV=maxVertex;

}


void LdpInstance::readGraph(std::ifstream& data,size_t maxVertex,char delim){
	std::string line;
	//	char delim = ' ';
	size_t lineCounter=0;
	std::getline(data, line);
	lineCounter++;
	std::cout << "called read graph" << std::endl;
	parameters.infoFile()<<"called read graph" << std::endl;
	std::vector<std::string> strings = split(line, delim);
	size_t numberOfVertices;

	if (strings.size() == 1) {
		if(parameters.isRestrictFrames()){
			numberOfVertices=maxVertex+3;
		}
		else{
			numberOfVertices = stoul(strings[0]);
			numberOfVertices += 2;  //to include s and t in the end of the list
		}
		s_ = numberOfVertices - 2;
		t_ = numberOfVertices - 1;

	} else {
		std::string str="first row must contain 1 number, detected ";
		str+=std::to_string(strings.size());
		str+="numbers";
		throw std::runtime_error(str);
	}



	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
	graph_ = andres::graph::Digraph<>(numberOfVertices);
	std::vector<double> inputCosts(numberOfVertices-2,parameters.getInputCost());
	std::vector<double> outputCosts(numberOfVertices-2,parameters.getOutputCost());


	// std::vector<std::pair<size_t,siz  Data<>e_t> > liftedEdges;
	vertexScore = std::vector<double>(numberOfVertices, 0);

	std::cout<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile()<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile().flush();
	//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 2) {
			throw std::runtime_error(
					std::string("Vertex and its score expected"));
		}


		unsigned int v = std::stoul(strings[0]);
		if(v>graph_.numberOfVertices()-3) continue;
		double score = std::stod(strings[1]);
		vertexScore[v] = score;

		if(strings.size()==4){
			inputCosts[v]=std::stod(strings[2]);
			outputCosts[v]=std::stod(strings[3]);
		}

	}

	for (int v = 0; v < numberOfVertices-2; ++v) {
		graph_.insertEdge(s_,v);
		edgeScore.push_back(inputCosts[v]);
		graph_.insertEdge(v,t_);
		edgeScore.push_back(outputCosts[v]);
	}

	size_t maxGap=parameters.getMaxTimeGapComplete();

	std::cout<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile().flush();
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 3) {
			throw std::runtime_error(
					std::string("Edge vertices and score expected, line "+std::to_string(lineCounter)));
		}

		unsigned int v = std::stoul(strings[0]);
		unsigned int w = std::stoul(strings[1]);

		if(v>numberOfVertices-3||w>numberOfVertices-3) continue;

		size_t gv=vertexGroups.getGroupIndex(v);
		size_t gw=vertexGroups.getGroupIndex(w);
		if(gw-gv>maxGap) continue;

		//if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
		double score = std::stod(strings[2]);
		auto edgeTest=graph_.findEdge(v,w);

		if(!edgeTest.first){  //if the edge does not exist
			graph_.insertEdge(v, w);
			edgeScore.push_back(score);

		}
		else{  //if the edge already exists, only update the score
			edgeScore[edgeTest.second]=score;

		}
	}
}

/*
template<typename EDGE_LABEL_ITERATOR>
	bool LdpInstance::check_feasiblity(EDGE_LABEL_ITERATOR begin, EDGE_LABEL_ITERATOR end) const{

		std::vector<bool> isOnPath(graph_.numberOfVertices(),0);
		std::vector<bool> isOnAnyPath(graph_.numberOfVertices(),0);

		//std::vector<size_t> pathToEnd(graph.numberOfVertices());



		size_t trackCounter=0;
		size_t activeEdgeCounter=0;
		size_t indexCounter=0;

		std::vector<double> solution0(graph_.numberOfEdges()+graphLifted_.numberOfEdges()+graphLifted_.numberOfVertices());
		for (auto it=begin;it!=end&&indexCounter<solution0.size();it++) {
				double value=*it;
				solution0[indexCounter]=value;
				indexCounter++;
		}

		for (int i = 0; i < graph_.numberOfEdgesToVertex(t_); ++i) {
			//std::cout<<i<<"-th edge to t"<<std::endl;

			size_t e=graph_.edgeToVertex(t_,i);
			size_t eVarIndex=getEdgeVarIndex(e);
			if(solution0[eVarIndex] > 0.5){
				//std::cout<<"check track "<<trackCounter<<std::endl;
				trackCounter++;
				activeEdgeCounter++;
				//std::map<size_t,size_t> vertexToDepth;
				size_t vertex=graph_.vertexToVertex(t_,i);
				isOnPath[vertex]=1;
				isOnPath[t_]=1;  //Maybe not necessary
				isOnAnyPath[vertex]=1;
				isOnAnyPath[t_]=1;

				std::list<size_t> path;
				path.push_front(t_);
				path.push_front(vertex);

				while(vertex!=s_){

					for (int j = 0; j < graphLifted_.numberOfEdgesFromVertex(vertex); ++j) {
						size_t le=graphLifted_.edgeFromVertex(vertex,j);
						size_t vertex2=graphLifted_.vertexFromVertex(vertex,j);
						size_t leVarIndex=getLiftedEdgeVarIndex(le);
						//size_t vertex2VarIndex=data_.getVertexVarIndex(vertex2);
						if(isOnPath[vertex2]&&solution0[leVarIndex]<0.5) {
							parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
							parameters.infoFile()<<"Error: Inactive lifted edge on an active path"<<std::endl;
							parameters.infoFile().flush();
							std::cout<<"check track "<<trackCounter<<std::endl;
							throw std::runtime_error(std::string("Inactive lifted edge on an active path"));
						}
						else if(!isOnPath[vertex2]&&solution0[leVarIndex]>0.5){
							parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
							parameters.infoFile()<<"Active lifted edge outside an active path"<<std::endl;
							parameters.infoFile().flush();
							std::cout<<"check track "<<trackCounter<<std::endl;
							throw std::runtime_error(std::string("Active lifted edge outside an active path"));
						}


					}


					size_t newVertex=0;
					bool nextFound=false;
					for (int j = 0; j < graph_.numberOfEdgesToVertex(vertex);
							++j) {
						e = graph_.edgeToVertex(vertex, j);
						eVarIndex = getEdgeVarIndex(e);

						if (solution0[eVarIndex] > 0.5) {
							newVertex = graph_.vertexToVertex(vertex, j);
							path.push_front(newVertex);
							isOnPath[newVertex]=1;
							isOnAnyPath[newVertex]=1;
							activeEdgeCounter++;
							if(nextFound){
								parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
								parameters.infoFile()<<"Multiple active incoming edges to vertex "<<vertex<<std::endl;
								parameters.infoFile().flush();
								std::cout<<"check track "<<trackCounter<<std::endl;
								std::cout<<"Multiple active incoming edges to vertex "<<vertex<<std::endl;
							}
							nextFound=true;


							//break;
						}
					}
					vertex=newVertex;



				}

				//std::cout<<std::endl;
				while(!path.empty()){
					isOnPath[path.front()]=0;
					path.pop_front();
				}
			}
		}

		for (int i = 0; i < graph_.numberOfVertices(); ++i) {
			if(i==s_||i==t_) continue;
			if((solution0[i]>0.5&&!isOnAnyPath[i])||(solution0[i]<0.5&&isOnAnyPath[i])){
				parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
				parameters.infoFile()<<"vertex "<<i<<" labeled "<<solution0[i]<<" is on Path "<<isOnAnyPath[i]<<std::endl;
				parameters.infoFile().flush();
				std::cout<<"check track "<<trackCounter<<std::endl;
				std::cout<<"vertex "<<i<<" labeled "<<solution0[i]<<" is on Path "<<isOnAnyPath[i]<<std::endl;
			}
		}
		size_t activeEdgeCounter2=0;
		for (int i = 0; i < graph_.numberOfEdges(); ++i) {
			size_t index=getEdgeVarIndex(i);
			size_t v=getVertexVarIndex(graph_.vertexOfEdge(i,0));
			size_t w=getVertexVarIndex(graph_.vertexOfEdge(i,1));

			if(v==s_||w==t_) continue;
			if(solution0[index]>0.5){
				activeEdgeCounter2++;
				//if(v==s||w==t) continue;
				if(solution0[v]<0.5||solution0[w]<0.5){
					parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
					parameters.infoFile()<<"edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
					parameters.infoFile().flush();
					std::cout<<"check track "<<trackCounter<<std::endl;
					std::cout<<"edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
				}
			}
			//		if(solution0[index]==0&&(solution0[v]==1||solution0[w]==1)){
			//			std::cout<<"edge "<<v<<","<<w<<" has label , but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
			//		}
		}
		for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
			size_t index=getLiftedEdgeVarIndex(i);
			size_t v=getVertexVarIndex(graphLifted_.vertexOfEdge(i,0));
			size_t w=getVertexVarIndex(graphLifted_.vertexOfEdge(i,1));
			if(v==s_||w==t_) continue;
			if(solution0[index]>0.5&&(solution0[v]<0.5||solution0[w]<0.5)){
				parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
				parameters.infoFile()<<"lifted edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
				parameters.infoFile().flush();
				std::cout<<"check track "<<trackCounter<<std::endl;
				std::cout<<"lifted edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
			}

		}

			std::cout<<"Solution checked. Tracks: "<<trackCounter<<std::endl;
			parameters.infoFile()<<"Solution checked. Tracks: "<<trackCounter<<std::endl;
			parameters.infoFile().flush();


}

*/







void LdpInstance::sparsifyBaseGraph(){
	std::cout<<"Sparsify base graph"<<std::endl;
	parameters.infoFile()<<"Sparsify base graph"<<std::endl;
	andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());
	std::vector<double> newBaseCosts;
	//std::vector<size_t> inOutEdges;
	size_t k=parameters.getKnnK();
	//std::vector<size_t> goodLongEdges;


	std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
	for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
		std::unordered_map<int,std::list<size_t>> edgesToKeep;
		size_t l0=vertexGroups.getGroupIndex(v0);
		for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
			size_t e=graph_.edgeFromVertex(v0,ne);
			size_t v1=graph_.vertexFromVertex(v0,ne);
			//			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
			if(v0==s_||v1==t_){

				finalEdges[e]=true;
			}
			else{
				size_t l1=vertexGroups.getGroupIndex(v1);
				size_t gap=l1-l0;
				if(gap<=parameters.getMaxTimeBase()){

					double cost=edgeScore[e];
					if(edgesToKeep.count(gap)>0){
						std::list<size_t>& smallList=edgesToKeep[gap];
						auto it=smallList.begin();
						double bsf=edgeScore[*it];
						//std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
						while(bsf>cost&&it!=smallList.end()){
							it++;
							size_t index=*it;
							if(it!=smallList.end()){
								bsf=edgeScore[index];
								//	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
							}
						}
						if(it!=smallList.begin()){
							smallList.insert(it,e);
							if(smallList.size()>k) smallList.pop_front();
						}
						else if(smallList.size()<k){
							smallList.push_front(e);
						}
					}
					else{
						edgesToKeep[gap].push_front(e);
					}
				}

			}
		}

		double bsf=0;
		for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
			if(edgesToKeep.count(gap)>0){
				auto& smallList=edgesToKeep[gap];
				for(size_t e:smallList){
					finalEdges[e]=true;
					if(edgeScore[e]<bsf){
						bsf=edgeScore[e];
					}
				}
			}
		}

		bool onlyImproving=parameters.isRequireImproving();
		for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
			if(edgesToKeep.count(gap)>0){
				auto& smallList=edgesToKeep[gap];
				for(size_t e:smallList){
					double score=edgeScore[e];
					if(score<=parameters.getBaseUpperThreshold()){
						finalEdges[e]=true;
					}
				}

			}
		}
	}

	for (int e = 0; e < graph_.numberOfEdges(); ++e) {
		if(finalEdges[e]){
			size_t v0=graph_.vertexOfEdge(e,0);
			size_t v1=graph_.vertexOfEdge(e,1);
			tempGraph.insertEdge(v0,v1);
			newBaseCosts.push_back(edgeScore[e]);
		}
	}

	if(newBaseCosts.size()!=tempGraph.numberOfEdges()){
		throw std::runtime_error("Error in base graph sparsification.");
	}



	graph_=tempGraph;
	edgeScore=newBaseCosts;

	if(graph_.numberOfEdges()!=newBaseCosts.size()){
		parameters.infoFile()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
		parameters.infoFile().flush();
		std::cout<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
	}
	else{
		std::cout<<"edge number and graph size match "<<std::endl;
		parameters.infoFile()<<"edge number and graph size match "<<std::endl;
		parameters.infoFile().flush();
	}


	if(parameters.getSmallIntervals()==0){
		//desc=initReachable(graph_,parameters,&vertexGroups);
		reachable=initReachableSet(graph_,&vertexGroups);
	}
	//desc=initReachable(graph_,parameters);
	reachable=initReachableSet(graph_);


	std::cout<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
	parameters.infoFile()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
	parameters.infoFile().flush();

}

void LdpInstance::sparsifyLiftedGraph(){

	std::cout<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile()<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile().flush();
	//TODO run automaticLifted to find candidates first

	double negMaxValue=0;
	double posMinValue=0;
	bool useAdaptive=false;
	//TODO adaptive lifted threshold in config file
	std::vector<double> newLiftedCosts;

	andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


	negMaxValue=parameters.getNegativeThresholdLifted();
	posMinValue=parameters.getPositiveThresholdLifted();


	std::unordered_map<size_t,std::set<size_t>> liftedEdges;
	for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
		std::unordered_set<size_t> alternativePath;
		for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
			size_t w=graph_.vertexFromVertex(v,i);
			for(size_t u:reachable[w]){
				if(u!=w) alternativePath.insert(u);
			}
		}
		for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
			size_t w=graphLifted_.vertexFromVertex(v,i);
			if(w!=t_){
				if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

			}
		}
	}


	std::cout<<"done"<<std::endl;
	parameters.infoFile()<<"done"<<std::endl;
	parameters.infoFile().flush();


	for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
		size_t v0=graphLifted_.vertexOfEdge(i,0);
		size_t v1=graphLifted_.vertexOfEdge(i,1);
		int l0=vertexGroups.getGroupIndex(v0);
		int l1=vertexGroups.getGroupIndex(v1);
		double cost=getLiftedEdgeScore(i);
		bool goodCost=(cost<negMaxValue)||(cost>posMinValue);
		if(isReachable(v0,v1)){

			int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
			bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
			if(timeConstraint&&goodCost){
				if(liftedEdges[v0].count(v1)>0){
					tempGraphLifted.insertEdge(v0,v1);
					newLiftedCosts.push_back(cost);
				}
				else{
					auto edgeTest=graph_.findEdge(v0,v1);
					if(edgeTest.first){
						edgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
					}

				}

			}
		}

	}



	liftedEdgeScore=newLiftedCosts;

	graphLifted_=tempGraphLifted;
	std::cout<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
	parameters.infoFile()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
	parameters.infoFile().flush();

	if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
		std::cout<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
		parameters.infoFile()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
	}
	else{
		std::cout<<"lifted edge number and lifted graph size match "<<std::endl;
		parameters.infoFile()<<"lifted edge number and lifted graph size match "<<std::endl;

	}
	parameters.infoFile().flush();

}



}}//End of namespaces
