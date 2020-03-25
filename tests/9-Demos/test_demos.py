import windse_driver
import glob, os, sys
import traceback

demo_path = "demo/documented/"

### Get all yaml files in the documented folder ###
dirList = glob.glob(demo_path+"**/*.yaml",recursive=True)

### list all strings that need to be filtered ###
filters = ["output","driver"]

### split path and filename ###
demos = []
for d in dirList:
    
    process_demo = True
    ### Filter out all yaml files in output/ folders ###
    for f in filters:
        process_demo = process_demo and (f not in d)

    if process_demo:
        path, name = os.path.split(d)
        demos.append([path,name])

        # print(demo_path[-1], demo_name[-1])

home_path = os.getcwd()

num_demos = len(demos)

### This allows specific demos to be run#
if len(sys.argv) > 1:
    test_list = [int(i) for i in sys.argv[1:]]
else:
    test_list = range(1,num_demos+1)
num_tests = len(test_list)

success_counter = 0
failed_test = []
test_results = []
for i, (demo_path, demo_name) in enumerate(demos):

    if i+1 in test_list:
        header_string = "========== Testing {:2d} of {:2d}: {:} ==========".format(i+1,num_demos,demo_name)
        sys.stdout.write("\x1b]2;"+header_string+"\x07")
        print("========================= Testing {:2d} of {:2d}: =========================".format(i+1,num_demos))
        print()
        print("Demo:     "+demo_path+"/"+demo_name)
        ### Switch to the demo directory ###
        os.chdir(demo_path)

        ### Run the test and save the results ###
        # save_stdout = sys.stdout
        # sys.stdout = open('trash', 'w')
        test_results.append(windse_driver.driver.test_demo(demo_name))
        # sys.stdout = save_stdout

        ### Return to the the original path ###
        os.chdir(home_path)
        print()
        print("===========================================================================")
        print()
        print()
        print()
        print()
        print()

for i, status in enumerate(test_results):

    (demo_path, demo_name) = demos[i]

    print("============================== Test {:2d} of {:2d}: =============================".format(i+1,num_demos))
    print()
    print("Demo:     "+demo_path+"/"+demo_name)

    if status[0]:
        print("Status:   Passed")
        print("Run Time: {: .2f} s".format(status[1]))
        success_counter += 1
    else:
        print("Status:   Failed")
        print("Reason:   "+repr(status[1].args[0]))
        print()
        print("-------------------------- Traceback ---------------------------")
        traceback.print_tb(status[2])
        print("-----------------------------------------------------------------")

        failed_test.append(i+1)    

    print()
    print("===========================================================================")
    print()
    print()
    print()

print("Testing Finished: {:d} of {:d} Passed".format(success_counter,num_tests))
if failed_test:
    print("Failing Tests: "+repr(failed_test))
