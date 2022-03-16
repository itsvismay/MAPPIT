import json
# files = ["TimevsIncreasingSizedAgentCircles/2_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/5_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/10_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/15_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/20_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/25_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/30_Agents/run0/agents.json",
# 		"TimevsIncreasingSizedAgentCircles/35_Agents/run0/agents.json"]

# files = ["TimevsIncreasingDofs8AgentCircle/run30/agents.json",
# 		"TimevsIncreasingDofs8AgentCircle/run60/agents.json",
# 		"TimevsIncreasingDofs8AgentCircle/run90/agents.json",
# 		"TimevsIncreasingDofs8AgentCircle/run120/agents.json",
# 		"TimevsIncreasingDofs8AgentCircle/run150/agents.json",
# 		"TimevsIncreasingDofs8AgentCircle/run180/agents.json"]

# files = ["./IncreasingEnvMeshRes/agents0.json",
# 		"IncreasingEnvMeshRes/agents1.json",
# 		"IncreasingEnvMeshRes/agents2.json",
# 		"IncreasingEnvMeshRes/agents3.json"]

# files = [r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\scaling_tests\2_agents\run0\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\scaling_tests\8_agents\run0\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\scaling_tests\10_agents\run0\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\scaling_tests\20_agents\run0\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\scaling_tests\30_agents\run0\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\airplane\3agents\run17\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\no_collisions\run-1\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\symmetric_collisions\run2\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\mass_asymmetric_collisions\run2\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\size_asymmetric_collisions\run7\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\size_mass_asymmetric_collisions\run9\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\three_agents\friends_asymmetric_grouping\run2\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\roomba_maze\scene_3\run2\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\bottleneck\bottleneck\run3\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\complex_maze\square_maze\five_agents\run11\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\circle_maze\seven_agents\run1\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\tunnel_maze\scene_1\run8\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\battlfield\agents.json',
# r'C:\Users\vismay\OneDrive - University of Toronto\CrowdsPngs\Timings\ricky_baboon_elephant\run-1\agents.json'
# ]
files = ["/Users/vismay/recode/crowds/Scenes/2_output_results/love/sparse/run28/agents.json"]

for i in range(0,len(files)):		
	data = json.load(open(files[i]))

	numAgents = len(data["agents"])
	numSegments = int(data["agents"][0]["segments"])
	preprocessTotalTime = float(data["timings"]["preprocess"])
	iterations = float(len(data["timings"]["iterations"]))
	avgTimePerIterationForSolve = [float(it["totalTime"])/int(it["output"]["funcCount"]) for it in data["timings"]["iterations"]] if isinstance(data["timings"]["iterations"], list) else [float(data["timings"]["iterations"]["totalTime"])/int(data["timings"]["iterations"]["output"]["funcCount"])]
	totalSolveTime = [float(it["totalTime"]) for it in data["timings"]["iterations"]] if isinstance(data["timings"]["iterations"], list) else [float(data["timings"]["iterations"]["totalTime"])]
	# print("----"+str(numAgents)+"------")
	print(numAgents, numAgents*numSegments, preprocessTotalTime/numAgents,sum(avgTimePerIterationForSolve), sum(totalSolveTime), len(totalSolveTime))
