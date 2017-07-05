import time
def ConnectivityMapper(fluid_nodes, solid_nodes):
    fluid_nodes=list(fluid_nodes)
    solid_nodes=list(solid_nodes)

    if len(fluid_nodes)!=len(solid_nodes):
        raise('Different length of interfaces')
    #~ print('Hola, estoy dentro')
    connectivity_list=[]
    id_list=[]
    tol=1e-8 #Geometric tolerance.
    i=0

    for node in fluid_nodes:
        x=node.X
        y=node.Y
        suptolx=x+tol
        inftolx=x-tol
        suptoly=y+tol
        inftoly=y-tol
        #~ print('Para nodo con coordenadas ',x,' ; ',y)
        #~ print(suptolx,' ',inftolx,' ',suptoly,' ',inftoly)
        j=0
        for node_solid in solid_nodes:
            #~ print(node_solid.Id)
            #~ err
            x_solid=node_solid.X
            y_solid=node_solid.Y

            if x_solid<suptolx and x_solid>inftolx and y_solid<suptoly and y_solid>inftoly:
                #~ print('Tenemos un caso!')
                connectivity_list.append([i,j])             # List w/ positions: i position in fluid corresponds to j position in solid
                id_list.append([node.Id,node_solid.Id])     # List w/ nodal id's: same concept but with nodal Id's.
                print('Fluid node: ',node.Id,'; Solid node: ',node_solid.Id)
                #~ time.sleep(3)
                break

            j+=1

        i+=1
    #~ print(connectivity_list)
    #~ err
    #~ time.sleep(15)
    return connectivity_list
