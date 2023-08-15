import numpy as np


class LEFTranslocator():
    
    def __init__(self, birthArray, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  numLEF):
        birthArray[0] = 0
        birthArray[len(birthArray)-1] = 0
        
        birthArray[stallProbLeft > 0.9] = 0        
        birthArray[stallProbRight > 0.9] = 0        
        
        self.numSite = len(birthArray)
        self.numLEF = numLEF
        self.stallLeft = stallProbLeft
        self.stallRight = stallProbRight
        self.falloff = deathProb
        self.pause = pauseProb
        self.birthProb = np.cumsum(birthArray, dtype=np.double)
        self.birthProb /= self.birthProb[-1]
        self.LEFs1 = np.zeros((self.numLEF), int)
        self.LEFs2 = np.zeros((self.numLEF), int)
        self.stalled1 = np.zeros(self.numLEF, int)
        self.stalled2 = np.zeros(self.numLEF, int)
        self.occupied = np.zeros(self.numSite, int)
        self.stallFalloff = stallFalloffProb
        self.occupied[0] = 1
        self.occupied[self.numSite - 1] = 1

        for ind in range(self.numLEF):
            self.birth(ind)


    def birth(self, ind):
    
        while True:
        
            pos = self.getStartSite()
            
            if pos >= self.numSite - 1:
                print("bad value", pos, self.birthProb[len(self.birthProb)-1])
                continue
                
            if pos <= 0:
                print("bad value", pos, self.birthProb[0])
                continue 
 
            if self.occupied[pos] == 1:
                continue
            
            self.LEFs1[ind] = pos
            self.LEFs2[ind] = pos
            
            self.occupied[pos] = 1
            
            if (pos < (self.numSite - 3)) and (self.occupied[pos+1] == 0):
                if np.random.random() > 0.5:
                    self.LEFs2[ind] = pos + 1
                    self.occupied[pos+1] = 1
            
            return


    def death(self):
    
        for i in range(self.numLEF):
        
            if self.stalled1[i] == 0:
                falloff1 = self.falloff[self.LEFs1[i]]
            else: 
                falloff1 = self.stallFalloff[self.LEFs1[i]]
                
            if self.stalled2[i] == 0:
                falloff2 = self.falloff[self.LEFs2[i]]
            else:
                falloff2 = self.stallFalloff[self.LEFs2[i]]              
            
            falloff = max(falloff1, falloff2)
            
            if np.random.random() < falloff:
                self.occupied[self.LEFs1[i]] = 0
                self.occupied[self.LEFs2[i]] = 0
                
                self.stalled1[i] = 0
                self.stalled2[i] = 0
                
                self.birth(i)
    
    
    def getStartSite(self):
    
        startSite = np.searchsorted(self.birthProb, np.random.random())
        
        return startSite
        

    def step(self):
        for i in range(self.numLEF):            
            stall1 = self.stallLeft[self.LEFs1[i]]
            stall2 = self.stallRight[self.LEFs2[i]]
                                    
            if np.random.random() < stall1:
                self.stalled1[i] = 1
                
            if np.random.random() < stall2:
                self.stalled2[i] = 1
                         
            cur1 = self.LEFs1[i]
            cur2 = self.LEFs2[i]
            
            if self.stalled1[i] == 0: 
                if self.occupied[cur1-1] == 0:
                    pause1 = self.pause[self.LEFs1[i]]
                    
                    if np.random.random() > pause1:
                        self.occupied[cur1 - 1] = 1
                        self.occupied[cur1] = 0
                        
                        self.LEFs1[i] = cur1 - 1
                        
            if self.stalled2[i] == 0:                
                if self.occupied[cur2 + 1] == 0:                    
                    pause2 = self.pause[self.LEFs2[i]]
                    
                    if np.random.random() > pause2:
                        self.occupied[cur2 + 1] = 1
                        self.occupied[cur2] = 0
                        
                        self.LEFs2[i] = cur2 + 1
        
        
    def steps(self,N):
        for i in range(N):
            self.death()
            self.step()
            
    def getOccupied(self):
        return np.array(self.occupied)
    
    def getLEFs(self):
        return np.array(self.LEFs1), np.array(self.LEFs2)
        
    def updateMap(self, cmap):
        cmap[self.LEFs1, self.LEFs2] += 1
        cmap[self.LEFs2, self.LEFs1] += 1

    def updatePos(self, pos, ind):
        pos[ind, self.LEFs1] = 1
        pos[ind, self.LEFs2] = 1

        
class LEFTranslocatorDynamicBoundary(LEFTranslocator):
     
    def __init__(self, birthArray, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  ctcfLifetime, ctcfUnboundTime, numLEF):
        self.stallLeft_init = np.copy(stallProbLeft)
        self.stallRight_init = np.copy(stallProbRight)        
        self.ctcfUnboundingRate = 1./ctcfLifetime if ctcfLifetime>0 else 0
        self.ctcfBoundingRate = 1./ctcfUnboundTime if ctcfUnboundTime>0 else 0

        super().__init__(birthArray, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  numLEF)

    def death_ctcf(self):     
        for i in range(self.numSite):
            if self.stallLeft[i] != 0:
                if np.random.random() < self.ctcfUnboundingRate:
                    self.stallLeft[i] = 0.0
                    for j in range(self.numLEF):
                        if i==self.LEFs1[j]:
                            self.stalled1[j] = 0
                    
            if self.stallRight[i]!=0:
                if np.random.random() < self.ctcfUnboundingRate:
                    self.stallRight[i]=0
                    for j in range(self.numLEF):
                        if i==self.LEFs2[j]:
                            self.stalled2[j] = 0
                          

    def birth_ctcf(self):        
        for i in range(self.numSite):
            if self.stallLeft[i] != self.stallLeft_init[i]:
                if np.random.random() < self.ctcfBoundingRate:
                    self.stallLeft[i] = self.stallLeft_init[i]
            if self.stallRight[i] != self.stallRight_init[i]:
                if np.random.random() < self.ctcfBoundingRate:
                    self.stallRight[i] = self.stallRight_init[i]
                        
                                
    def steps(self, N):
        for i in range(N):
            self.death_ctcf()
            self.birth_ctcf()
            
            self.death()
            self.step()
            
        
        