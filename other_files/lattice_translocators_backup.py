import numpy as np


class LEFTranslocator():
    
    def __init__(self, emissionProb, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  numLEF):
        emissionProb[0] = 0
        emissionProb[len(emissionProb)-1] = 0
        
        emissionProb[stallProbLeft > 0.9] = 0        
        emissionProb[stallProbRight > 0.9] = 0        
        
        self.N = len(emissionProb)
        self.M = numLEF
        self.emission = emissionProb
        self.stallLeft = stallProbLeft
        self.stallRight = stallProbRight
        self.falloff = deathProb
        self.pause = pauseProb
        cumem = np.cumsum(emissionProb, dtype=np.double)
        cumem /= cumem[-1]
        self.cumEmission = cumem
        self.LEFs1 = np.zeros((self.M), int)
        self.LEFs2 = np.zeros((self.M), int)
        self.stalled1 = np.zeros(self.M, int)
        self.stalled2 = np.zeros(self.M, int)
        self.occupied = np.zeros(self.N, int)
        self.stallFalloff = stallFalloffProb
        self.occupied[0] = 1
        self.occupied[self.N - 1] = 1
        self.maxss = 1000000
        self.curss = 99999999

        for ind in range(self.M):
            self.birth(ind)


    def birth(self, ind):
    
        while True:
        
            pos = self.getss()
            
            if pos >= self.N - 1:
                print("bad value", pos, self.cumEmission[len(self.cumEmission)-1])
                continue
                
            if pos <= 0:
                print("bad value", pos, self.cumEmission[0])
                continue 
 
            if self.occupied[pos] == 1:
                continue
            
            self.LEFs1[ind] = pos
            self.LEFs2[ind] = pos
            
            self.occupied[pos] = 1
            
            if (pos < (self.N - 3)) and (self.occupied[pos+1] == 0):
                if np.random.random() > 0.5:
                    self.LEFs2[ind] = pos + 1
                    self.occupied[pos+1] = 1
            
            return


    def death(self):
    
        for i in range(self.M):
        
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
    
    
    def getss(self):
    
        if self.curss >= self.maxss - 1:
            foundArray = np.array(np.searchsorted(self.cumEmission, np.random.random(self.maxss)), dtype = np.long)
            self.ssarray = foundArray

            self.curss = -1
        
        self.curss += 1
        
        return self.ssarray[self.curss]
        

    def step(self):
        for i in range(self.M):            
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
     
    def __init__(self, emissionProb, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  ctcfLifetime, ctcfUnboundTime, numLEF):
        self.stallLeft_init = np.copy(stallProbLeft)
        self.stallRight_init = np.copy(stallProbRight)        
        self.ctcfUnboundingRate = 1./ctcfLifetime if ctcfLifetime>0 else 0
        self.ctcfBoundingRate = 1./ctcfUnboundTime if ctcfUnboundTime>0 else 0

        super().__init__(emissionProb, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  numLEF)

    def death_ctcf(self):     
        for i in range(self.N):
            if self.stallLeft[i] != 0:
                if np.random.random() < self.ctcfUnboundingRate:
                    self.stallLeft[i] = 0.0
                    for j in range(self.M):
                        if i==self.LEFs1[j]:
                            self.stalled1[j] = 0
                    
            if self.stallRight[i]!=0:
                if np.random.random() < self.ctcfUnboundingRate:
                    self.stallRight[i]=0
                    for j in range(self.M):
                        if i==self.LEFs2[j]:
                            self.stalled2[j] = 0
                          

    def birth_ctcf(self):        
        for i in range(self.N):
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
            
        
        