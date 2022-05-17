from coffea import processor, hist
import awkward as ak


class CountEvents(processor.ProcessorABC):
    '''Coffea processor that accumulates the original sum of weights
    for different NanoAOD samples with pre-selections, which are stored 
    in the  "Runs" tree. 
    '''
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            'sumw':processor.defaultdict_accumulator(float)
        })

    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata['dataset']
        if dataset != 'Data':
            output['sumw'][dataset] += ak.sum(events.genEventSumw)
            # print(dataset)
            # print( ak.sum(events.genEventSumw))

            # output['sumw'][dataset+'_sing'] += events.genEventSumw[0]
            # output['sumw'][dataset+'_2highest'] += events.genEventSumw[0]
            # output['sumw'][dataset+'_+missingpT'] += events.genEventSumw[0]
        else:
            output['sumw'][dataset] += len(events)
            print(dataset)

            # output['sumw'][dataset+'_sing'] += len(events)
            # output['sumw'][dataset+'_2highest'] += len(events)
            # output['sumw'][dataset+'_+missingpT'] += len(events)

        return output
    
    def postprocess(self, accumulator):
        return accumulator
