from dolfin import Constant
# from windse import windse_parameters

# current_tab = windse_parameters.current_tab

def BaseHeight(x,y,ground,dx=0,dy=0):
    return Constant(ground(float(x),float(y),dx=dx,dy=dx))

# def fprint(self,string,tab=None,offset=0,special=None):
#     """
#     This is just a fancy print function that will tab according to where
#     we are in the solve

#     Args:
#         string (str): the string for printing

#     :Keyword Arguments:
#         * **tab** (*int*): the tab level

#     """
#     ### Check Processor ###
#     rank = 0
#     if rank == 0:
#         ### Check if tab length has been overridden
#         if tab is None:
#             tab = self.current_tab
        
#         ### Check if we are starting or ending a section
#         if special=="header":
#             self.current_tab += 1
#             self.fprint("",tab=tab)
#         elif special =="footer":
#             self.current_tab -= 1
#             tab -= 1
#             self.fprint("",tab=tab+1)

#         ### Apply Offset if provided ###
#         tab += offset

#         ### Create Tabbed string
#         tabbed = "|    "*tab

#         ### Apply Tabbed string
#         if isinstance(string,str):
#             string = tabbed+string
#         else:
#             string = tabbed+repr(string)

#         ### Print
#         print(string)
#         # sys.stdout.flush()

#         if special=="header":
#             self.fprint("",tab=tab+1)
