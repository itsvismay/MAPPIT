import json 

files = ["TimevsIncreasingSizedAgentCircles/2_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/5_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/10_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/15_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/20_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/25_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/30_Agents/run0/agents.json",
		"TimevsIncreasingSizedAgentCircles/35_Agents/run0/agents.json"]

files = ["TimevsIncreasingDofs8AgentCircle/run30/agents.json",
		"TimevsIncreasingDofs8AgentCircle/run60/agents.json",
		"TimevsIncreasingDofs8AgentCircle/run90/agents.json",
		"TimevsIncreasingDofs8AgentCircle/run120/agents.json",
		"TimevsIncreasingDofs8AgentCircle/run150/agents.json",
		"TimevsIncreasingDofs8AgentCircle/run180/agents.json"]

files = ["IncreasingEnvMeshRes/agents0.json",
		"IncreasingEnvMeshRes/agents1.json",
		"IncreasingEnvMeshRes/agents2.json",
		"IncreasingEnvMeshRes/agents3.json"]

for i in range(0,len(files)):		
	data = json.load(open(files[i]))

	numAgents = len(data["agents"])
	numSegments = data["agents"][0]["segments"]
	preprocessTotalTime = data["timings"]["preprocess"]
	iterations = float(len(data["timings"]["iterations"]))
	envMeshBdrySize = float(data["env"]["sizeV"][0])
	avgTimePerIterationForSolve = [float(it["totalTime"])/int(it["output"]["funcCount"]) for it in data["timings"]["iterations"]]
	# print("----"+str(numAgents)+"------")
	#print(numAgents)
	# print(numSegments)
	# print(preprocessTotalTime)
	# print(iterations)
	print(sum(avgTimePerIterationForSolve))
	# print("----------")