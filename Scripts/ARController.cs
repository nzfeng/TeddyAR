using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using GoogleARCore;
using GoogleARCoreInternal;

public class ARController : MonoBehaviour
{
    public GameObject ARCamera;
    public SketchController sketchController;
    public List<Vector3> points = new List<Vector3>();
	public bool isSketching;

    // Start is called before the first frame update
    void Start()
    {
		isSketching = true;

        sketchController = ARCamera.GetComponent<SketchController>();
    }

    // Update is called once per frame
    void Update()
    {

        // Check ARCore session status.
        if (Session.Status != SessionStatus.Tracking)
        {
            int lostTrackingSleepTimeout = 15;
            Screen.sleepTimeout = lostTrackingSleepTimeout;
            return;
        }
        Screen.sleepTimeout = SleepTimeout.NeverSleep;

        if (Input.touchCount < 1 && Input.GetTouch(0).phase != TouchPhase.Began)
        {
            Screen.sleepTimeout = SleepTimeout.NeverSleep;
            return;
        }
        else if (Input.touchCount > 0 && Input.GetTouch(0).phase == TouchPhase.Ended)
        {
            if (isSketching)
            {
                Vector3 camPos = ARCamera.transform.position;
                Vector3 camDirection = ARCamera.transform.forward;
                float spawnDistance = 0.7f;
                Vector3 spawnPos = camPos + (camDirection * spawnDistance);
                Debug.Log("spawnPos: " + spawnPos);
                sketchController.spawnPos = spawnPos;

                if (points.Count > 0)
                {
                    isSketching = false;
                    sketchController.isSketching = false;
                }
                points.Clear();
            }
        }
        else if (Input.touchCount > 0 && Input.GetTouch(0).phase == TouchPhase.Began)
        {
            isSketching = true;
            sketchController.isSketching = true;
            sketchController.points.Clear();
        }
        else if (Input.touchCount > 0 && Input.GetTouch(0).phase == TouchPhase.Moved)
        { 
            Touch touch = Input.touches[0];
            Vector3 tappoint = new Vector3((touch.position.x-Screen.width/2.0f)/Screen.width, (touch.position.y-Screen.height/2.0f)/Screen.height, 0.0f); // screen space
            // Debug.Log(touch.position.x / Screen.width + " " + touch.position.y / Screen.height); // confirm screen coordinates
            points.Add(tappoint);
            // Deep copy
            sketchController.points.Add(new Vector3(tappoint.x, tappoint.y, tappoint.z));
        }

        // Check if the user has touched any of the tracked planes.
        TrackableHit hit;
        TrackableHitFlags raycastFilter = TrackableHitFlags.PlaneWithinBounds | TrackableHitFlags.PlaneWithinPolygon;
        // if (Frame.Raycast(touch.position.x, touch.position.y, raycastFilter, out hit))
        // {
		// 	// Place the portal on top of the plane that we touched.
		// 	flag = true;
		// 	// Create new anchor.
		// 	Anchor anchor = hit.Trackable.CreateAnchor(hit.Pose);
		// 	// Set the position of the portal to be the same as the hit position.
		// 	Portal.transform.position = hit.Pose.position;
		// 	Portal.transform.rotation = hit.Pose.rotation;
		// 	// Have portal face the camera. Only rotate portal around the y-axis.
		// 	Vector3 cameraPosition = ARCamera.transform.position;
		// 	cameraPosition.y = hit.Pose.position.y;
		// 	Portal.transform.LookAt(cameraPosition, transform.up);
		// 	// Attach portal to anchor.
		// 	Portal.transform.SetParent (anchor.transform);
		// 	// Display.
		// 	Portal.SetActive(true);
        // }
        


    }
}
