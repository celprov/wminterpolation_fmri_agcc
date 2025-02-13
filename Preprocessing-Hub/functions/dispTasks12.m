function dispTasks(tasksTodo, tasksDone)
    global OLDSEGMENT;
    global SEGMENT;
    global SEGMENTDARTEL;

    disp('The following tasks have already been done:')
    if tasksDone.realign, disp('-Realignment');  end
    if isfield(tasksDone, 'QC') && tasksDone.QC, disp('-Quality Control'); end
    if tasksDone.coregister, disp('-Coregistration'); end
    if tasksDone.coregReslice, disp('-coregReslice'); end
    if tasksDone.smooth, disp('-Smooth'); end
    if tasksDone.segment == SEGMENTDARTEL, disp('-Segmentation + Dartel'); end
    if tasksDone.segment == SEGMENT, disp('-Segmentation'); end
    if tasksDone.segment == OLDSEGMENT, disp('-Old Segmentation'); end
    if tasksDone.label, disp('-Labeling'); end
    
    disp('The following tasks have to be done:')
    if tasksTodo.realign, disp('-Realignment');  end
    if tasksTodo.QC, disp('-Quality Control'); end
    if tasksTodo.coregister, disp('-Coregistration'); end
    if tasksTodo.coregReslice, disp('-coregReslice'); end
    if tasksTodo.smooth, disp('-Smooth'); end
    if tasksTodo.segment == SEGMENTDARTEL, disp('-Segmentation + Dartel'); end
    if tasksTodo.segment == SEGMENT, disp('-Segmentation'); end
    if tasksTodo.segment == OLDSEGMENT, disp('-Old Segmentation'); end
    if tasksTodo.label, disp('-Labeling'); end
return
    